#include "SatNavRel.hpp"
#include "SatNav.hpp"

SatNavRel::SatNavRel(SatNav& passive, SatNav& active) : pas(passive), act(active) {
    pas.solve(); // решать одновременно?

    // тру стейтс
}

void SatNavRel::solve_relative(char et, int ti, int tf) {

    error_type = et;

    for (const auto& raw_mg_pas : pas.raw_measurements_groupped) {

        if ((ti > 0 || tf > 0) && (raw_mg_pas.time < ti || raw_mg_pas.time > tf)) continue;

        RawMeasurementGroupped raw_mg_act = act.get_raw_measurement_groupped_at(raw_mg_pas.time);

        if (raw_mg_act.time == -1) continue;

        std::vector<RefinedMeasurement> ref_ms(32);

        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            unsigned prn_index = prn_id - 1;

            if (!pas.check_raw(raw_mg_pas.raw_measurements[prn_index]) || !act.check_raw(raw_mg_act.raw_measurements[prn_index])) continue;

            RefinedMeasurement ref_m = refine_raw(raw_mg_pas.raw_measurements[prn_index], raw_mg_act.raw_measurements[prn_index]);

            ref_ms[prn_index] = ref_m;
        }

        RefinedMeasurementGroupped ref_mg = {raw_mg_pas.time, ref_ms};
        refined_measurements_groupped.push_back(ref_mg);

        SolutionState solution = calculate_solution(ref_mg);

        // if (error_type != 'L') {
        //     if (solution.is_solved) {
        //         std::vector<unsigned> low = check_low(solution, ref_mg);
        //         if (low.size() > 0) {
        //             for (unsigned prn_id : low) {
        //                 unsigned prn_index = prn_id - 1;
        //                 ref_mg.refined_measurements[prn_index].is_present = false;
        //             }
        //             solution = calculate_solution(ref_mg);
        //         }
        //     }
        // }

        solution_states.push_back(solution);
    }
   
    error_type = '0';
}

RefinedMeasurement SatNavRel::refine_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act) {
    double pseudorange_delta = raw_m_pas.L1_range - raw_m_act.L1_range;
    // Вычисляем положение НКА очень грубо
    // Может вычесть хотя бы ошибку часов НКА?
    std::vector<double> gps_position = pas.handler.get_state(raw_m_pas.prn_id, raw_m_pas.time).position;

    return {true, raw_m_pas.time, raw_m_pas.prn_id, pseudorange_delta, 0, gps_position};
}

SolutionState SatNavRel::calculate_solution(const RefinedMeasurementGroupped& ref_mg) const {
    SolutionState solution;
    solution.time = ref_mg.time;

    std::vector<double> solution_pas = pas.get_solution_state_at(ref_mg.time).position;
    std::vector<double> dPR;
    std::vector<std::vector<double>> X;
    size_t n = 0;

    for (const auto& ref_m : ref_mg.refined_measurements) {
        if (ref_m.is_present) {
            dPR.push_back(ref_m.pseudorange);
            X.push_back(ref_m.gps_position);
            n++;
        }
    }

    std::vector<double> U(n);
    Matrix B(n, 3);

    for (size_t i = 0; i < n; i++) {
        size_t j = i < n - 1 ? i + 1 : 0;
        std::vector<double> DXi = X[i] - solution_pas;
        std::vector<double> DXj = X[j] - solution_pas;
        double Di = abs(DXi);
        double Dj = abs(DXj);

        //delta

        U[i] = dPR[i] - dPR[j];
        for (size_t k = 0; k < 3; k++) {
            B.at(i, k) = DXi[k] / Di -  DXj[k] / Dj;
        }
    }

    Matrix B1(3, 3);
    try {
        B1 = (B.transpose() * B).inverse();
    } catch (std::runtime_error re) {
        solution.is_solved = false;
        solution.failure_type = 'I';
        return solution;
    }

    // GDOP

    std::vector<double> dX = B1 * B.transpose() * U;

    std::cout << abs(dX) << std::endl;

    solution.position = dX;
    solution.is_solved = true;

    return solution;
}

const State& SatNavRel::get_true_state_at(int time) const {
    return *get_true_state_iterator(time);
}

const std::vector<SolutionState>& SatNavRel::get_solution_states() const {
    return solution_states;
}

std::vector<State>::const_iterator SatNavRel::get_true_state_iterator(int time) const {
    size_t n = true_states.size();

    size_t lo = 0, hi = n - 1;
    while (lo <= hi) {
        size_t mid = lo + (hi - lo) / 2;

        if (true_states[mid].time == time) {
            return true_states.cbegin() + mid;
        }

        if (true_states[mid].time < time) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    return true_states.cend();
}
