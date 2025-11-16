#include "SatNav.hpp"

SatNav::SatNav(const std::string& gnv_filename, const std::string& gps_filename, 
               const GPSHandler& handler) : handler(handler) {
    true_states = DataParser::load_grace_fo_gnv_data(gnv_filename);
    raw_measurements_groupped = DataParser::load_grace_fo_gps_data(gps_filename);
}

SatNav::SatNav(const SatNav& sn) : handler(sn.handler), 
                                   true_states(sn.true_states), 
                                   raw_measurements_groupped(sn.raw_measurements_groupped) { }

void SatNav::solve(char et, unsigned ti, unsigned tf) {     // Очень мало точек, надо разобраться
    logger.log("Beginning to solve...");

    error_type = et;

    for (auto raw_mg_it = raw_measurements_groupped.begin();
         raw_mg_it != raw_measurements_groupped.end();
         raw_mg_it++) {
        
        RawMeasurementGroupped raw_mg = *raw_mg_it;

        if ((ti > 0 || tf > 0) && (raw_mg.time < ti || raw_mg.time > tf)) {
            continue;
        }

        // logger.log("t = " + std::to_string(raw_mg.time));

        std::vector<RefinedMeasurement> ref_ms(32);

        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            unsigned prn_index = prn_id - 1;

            RawMeasurement raw_m = raw_mg.raw_measurements[prn_index];

            if (!raw_m.is_present) {
                continue;
            }

            if (raw_m.qualflg != 0) {
                continue;
            }

            if (error_type != 'S') {
                if (raw_m.L1_SNR < SNR_threshold || raw_m.L2_SNR < SNR_threshold) {
                    continue;
                }
            }
            
            if (error_type != 'F') {
                if (check_fading(raw_m, raw_mg_it)) {
                    continue;
                }
            }

            RefinedMeasurement ref_m = refine_raw(raw_m);

            if (error_type != 'H') {
                ref_m = hatch_filter(ref_m);
            }

            ref_ms[prn_index] = ref_m;
        }

        RefinedMeasurementGroupped ref_mg = {raw_mg.time, ref_ms};
        refined_measurements_groupped.push_back(ref_mg);

        SolutionState solution = calculate_solution(ref_mg);

        if (error_type != 'L') {
            if (solution.is_solved) {
                std::vector<unsigned> low = check_low(solution, ref_mg);
                if (low.size() > 0) {
                    for (unsigned prn_id : low) {
                        unsigned prn_index = prn_id - 1;
                        ref_mg.refined_measurements[prn_index].is_present = false;
                    }
                    solution = calculate_solution(ref_mg);
                }
            }
        }
        // if (!solution.is_solved) {
        //     logger.log("Failed due to " + std::string(1, solution.failure_type));
        //     if (solution.failure_type == 'G') logger.log("GDOP = " + std::to_string(solution.GDOP));
        // }
        solution_states.push_back(solution);
    }
   
    error_type = '0';

    logger.log("Finished solving");
}

RefinedMeasurement SatNav::refine_raw(const RawMeasurement& raw_m) {
    RefinedMeasurement ref_m1 = apply_clock_and_relativistic_errors(raw_m, 1);
    RefinedMeasurement ref_m2 = apply_clock_and_relativistic_errors(raw_m, 2);

    double pseudorange = ref_m1.pseudorange;
    std::vector<double> gps_position = ref_m1.gps_position;

    if (error_type != 's') {
        double phi = ref_m1.pseudorange / c * earth_rotation_rate;
        Matrix rot = rotation(-phi, 'z');
        gps_position = rot * gps_position;
    }

    if (error_type != 'I') {
        pseudorange = ref_m1.pseudorange * C1 + ref_m2.pseudorange * C2;
    }
    double carrier_phase = raw_m.L1_phase * C1 + raw_m.L2_phase * C2;

    return {true, raw_m.time, raw_m.prn_id, pseudorange, carrier_phase, gps_position};
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

    if (error_type != 'R') {
        double relativistic_error = gs.relativistic_error;
        delay += relativistic_error;
    }

    return {true, raw_m.time, raw_m.prn_id, delay * c, 0, gs.position};
}

SolutionState SatNav::calculate_solution(const RefinedMeasurementGroupped& ref_mg) const {
    SolutionState solution;
    solution.time = ref_mg.time;

    std::vector<double> PR;
    std::vector<std::vector<double>> X;
    unsigned n = 0;

    for (const auto& ref_m : ref_mg.refined_measurements) {
        if (ref_m.is_present) {
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

bool SatNav::check_fading(const RawMeasurement& raw_m,
                  const std::vector<RawMeasurementGroupped>::iterator& raw_mg_it) {
    unsigned t0 = raw_m.time;
    unsigned prn_index = raw_m.prn_id - 1;

    auto it_fwd = raw_mg_it;
    while(true) {
        it_fwd++;
        if (it_fwd == raw_measurements_groupped.end()) return true;

        unsigned t = it_fwd->time;
        if (t - t0 > fadeout_time) break;

        if (!it_fwd->raw_measurements[prn_index].is_present) return true;
    }

    auto it_bwd = raw_mg_it;
    while(true) {
        if (it_fwd == raw_measurements_groupped.begin()) return true;
        it_fwd--;

        unsigned t = it_fwd->time;
        if (t0 - t > fadeout_time) break;

        if (!it_fwd->raw_measurements[prn_index].is_present) return true;
    }
    
    return false;
}

std::vector<unsigned> SatNav::check_low(const SolutionState& solution, const RefinedMeasurementGroupped& ref_mg) {
    std::vector<unsigned> low;

    for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
        unsigned prn_index = prn_id - 1;
        RefinedMeasurement ref_m = ref_mg.refined_measurements[prn_index];

        if (!ref_m.is_present) continue;
        
        std::vector<double> gps_relative = ref_m.gps_position - solution.position;
        double zenith_angle = angle_between(solution.position, gps_relative) * 180.0 * M_1_PI;
        if (90.0 - zenith_angle < mask_angle) {
            low.push_back(prn_id);
        }
    }

    return low;
}

RefinedMeasurement SatNav::hatch_filter(const RefinedMeasurement& ref_m) {
    unsigned prn_index = ref_m.prn_id - 1;

    RefinedMeasurement ref_m_hatch = ref_m;

    if (refined_measurements_groupped.size() > 0) {

        RefinedMeasurementGroupped ref_mg_last = *(--refined_measurements_groupped.end());  

        if (ref_mg_last.refined_measurements[prn_index].is_present) {
            double pseudorange_prev = ref_mg_last.refined_measurements[prn_index].pseudorange;
            double carrier_phase_prev = ref_mg_last.refined_measurements[prn_index].carrier_phase;
            double delta_phase = ref_m.carrier_phase - carrier_phase_prev;

            ref_m_hatch.pseudorange = hatch_constant * ref_m.pseudorange +
                                      (1 - hatch_constant) * (pseudorange_prev + delta_phase);
        }
    }

    return ref_m_hatch;
}

const State& SatNav::get_true_state_at(unsigned time) {
    unsigned n = true_states.size();

    unsigned lo = 0, hi = n - 1;
    while (lo <= hi) {
        unsigned mid = lo + (hi - lo) / 2;

        if (true_states[mid].time == time) {
            return true_states[mid];
        }

        if (true_states[mid].time < time) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    return empty_state;
}

const SolutionState& SatNav::get_solution_state_at(unsigned time) const {
    unsigned n = solution_states.size();

    unsigned lo = 0, hi = n - 1;
    while (lo <= hi) {
        unsigned mid = lo + (hi - lo) / 2;
        if (solution_states[mid].time == time) {
            return solution_states[mid];
        }

        if (solution_states[mid].time < time) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    return empty_solution_state;
}

const std::vector<SolutionState>& SatNav::get_solution_states() const {
    return solution_states;
}

const std::vector<RefinedMeasurementGroupped>& SatNav::get_refined_measurements_groupped() const {
    return refined_measurements_groupped;
}
