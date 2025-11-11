#include "SatNav.hpp"

SatNav::SatNav(const std::string& gnv_filename, const std::string& gps_filename, 
               const GPSHandler& handler) : handler(handler) {
    
    true_states = DataParser::load_grace_fo_gnv_data(gnv_filename);
    raw_measurements_groupped = DataParser::load_grace_fo_gps_data(gps_filename);
}

// мб и не войд, надо подумать
void SatNav::solve(unsigned ti, unsigned tf, char et) {
    error_type = et;

    for (const auto& raw_mg : raw_measurements_groupped) {

        if ((ti > 0 || tf > 0) && (raw_mg.time < ti || raw_mg.time > tf)) {
            continue;
        }

        std::vector<RefinedMeasurement> ref_ms(32);

        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            unsigned prn_index = prn_id - 1;

            RawMeasurement raw_m = raw_mg.raw_measurements[prn_index];
            
            // fading
            // qualflg
            // snr

            RefinedMeasurement ref_m = refine_raw(raw_m);
            ref_ms[prn_index] = ref_m;

            // hatch
        }

        RefinedMeasurementGroupped ref_mg = {raw_mg.time, ref_ms};
        refined_measurements_groupped.push_back(ref_mg);

        SolutionState solution = calculate_solution(ref_mg);

        if (solution.is_solved) {

            // low

            solution_states.push_back(solution);
        } else {
            // solution.failure_type
        }
    }
   
    error_type = '0';
}

RefinedMeasurement SatNav::refine_raw(const RawMeasurement& raw_m) {
    RefinedMeasurement ref_m1 = apply_clock_and_relativistic_errors(raw_m, 1);
    RefinedMeasurement ref_m2 = apply_clock_and_relativistic_errors(raw_m, 2);

    double phi = ref_m1.pseudorange / c * earth_rotation_rate;
    Matrix rot = rotation(-phi, 'z');
    std::vector<double> gps_position = rot * ref_m1.gps_position;
    
    double pseudorange = ref_m1.pseudorange * C1 + ref_m2.pseudorange * C2;

    return {raw_m.time, raw_m.prn_id, pseudorange, gps_position};
}

RefinedMeasurement SatNav::apply_clock_and_relativistic_errors(const RawMeasurement& raw_m, unsigned frequency) {
    /*  Эффект                     t       t * v    t * c   
    (1) Задержка распространения - 80 мс - 240 км -  ---
    (2) Ошибка часов НКА         - 3 мс  - 10 м   - 1000 км
    (3) Релятивизм               - 20 нс - 60 мкм - 6 м

    При расчете ошибки часов учитываем только (1)
    При расчете эфемерид учитываем (1), (2)
    При корректировке псевдодальности учитываем (1), (2), (3) */
    
    double pseudorange;
    if (frequency == 1) {
        pseudorange = raw_m.L1_range;
    } else if (frequency == 2) {
        pseudorange = raw_m.L2_range;
    }

    double delay = 0;

    double propagation_delay = pseudorange / c;
    delay += propagation_delay;

    double clock_error = handler.get_clock_error(raw_m.prn_id, raw_m.time - delay);
    delay += clock_error;

    GPSState gs = handler.get_state(raw_m.prn_id, raw_m.time - delay);

    double relativistic_error = gs.relativistic_error;
    delay += relativistic_error;

    return {raw_m.time, raw_m.prn_id, delay * c, gs.position};
}

SolutionState SatNav::calculate_solution(const RefinedMeasurementGroupped& ref_mg) const {

    SolutionState solution;

    std::vector<double> PR;
    std::vector<std::vector<double>> X;
    unsigned n = 0;

    for (const auto& ref_m : ref_mg.refined_measurements) {
        if (ref_m.prn_id != 0) {
            PR.push_back(ref_m.pseudorange);
            X.push_back(ref_m.gps_position);
            n++;
        }
    }

    double eps = 0.5;
    std::vector<double> x = {0.0, 0.0, 0.0};
    double c_tau = 0.0;

    std::vector<double> U(n);
    Matrix B(n, 4, 1.0);

    while (true) {
        for (unsigned i = 0; i < n; i++) {
            std::vector<double> DX = X[i] - x;
            double D = abs(DX);

            U[i] = PR[i] - D;
            for (unsigned j = 0; j < 3; j++) {
                B.at(i, j) = DX[j] / D;
            }
        }
        
        Matrix B1(4, 4);
        try {
            B1 = (B.transpose() * B).inverse();
        } catch (std::runtime_error re) {
            solution.is_solved = false;
            solution.failure_type = 'I';
            return solution;
        }

        double GDOP = sqrt(B1.trace());
        solution.GDOP = GDOP;

        if (GDOP > GDOP0) {
            solution.is_solved = false;
            solution.failure_type = 'G';
            return solution;
        }

        std::vector<double> dX = B1 * B.transpose() * U * (-1.0);

        std::vector<double> dX_x = std::vector<double>(dX.begin(), dX.begin() + 3);
        double dX_c_tau = dX[3];
    
        double delta = abs(dX_x) + std::abs(dX_c_tau - c_tau);
        if (delta < eps) {
            break;
        }

        x = x + dX_x;
        c_tau = dX_c_tau;
    }

    solution.time = ref_mg.time;
    solution.position = x;
    solution.is_solved = true;

    return solution;
}

