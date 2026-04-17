#include "SatNav.hpp"
#include "SatNavRel.hpp"

// Собирать полный путь либо здесь, либо в DataParser

SatNav::SatNav(const Date& date, const GRACE_SATS sat, const GPSHandler& handler) : handler(handler) {
    bool is_fo = sat == GRACE_SATS::C || sat == GRACE_SATS::D;
    if (is_fo) {
        true_states = DataParser::load_grace_fo_gnv_data(date, sat);
        raw_measurements_groupped = DataParser::load_grace_fo_gps_data(date, sat);
    } else {
        true_states = DataParser::load_grace_gnv_data(date, sat);
        raw_measurements_groupped = DataParser::load_grace_gps_data(date, sat);
    }
}

SatNav::SatNav(const Date& date, const SWARM_SATS sat, const GPSHandler& handler) : handler(handler) {

    true_states = DataParser::load_swarm_nav_data(date, sat);
    raw_measurements_groupped = DataParser::load_swarm_gps_data(date, sat);
}

SatNav::SatNav(const SatNav& sn) : handler(sn.handler), 
                                   true_states(sn.true_states), 
                                   raw_measurements_groupped(sn.raw_measurements_groupped),
                                   acceleration_measurements(sn.acceleration_measurements) {}

void SatNav::solve(char et, unsigned ti, unsigned tf) {
    logger.log("Beginning to solve...");

    error_type = et;

    unsigned t0 = raw_measurements_groupped[0].time;

    for (const auto& raw_mg : raw_measurements_groupped) {
        
        if (ti > 0 || tf > 0) {
            if ((raw_mg.time - t0) < ti || (raw_mg.time - t0) >= tf) continue;
        }

        logger.logv("Current time", raw_mg.time);

        std::vector<RefinedMeasurement> ref_ms(32);
        for (unsigned prn_id = 1; prn_id <= 32; ++prn_id) {
            unsigned prn_index = prn_id - 1;

            if (!check_raw(raw_mg.raw_measurements[prn_index])) continue;

            logger.logv("Refining", prn_id);

            RefinedMeasurement ref_m = refine_raw(raw_mg.raw_measurements[prn_index]);

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

        logger.log("Is solved: " + std::to_string(solution.is_solved) + ' ' + solution.failure_type);

        if (solution.is_solved) {
            auto ts = get_true_state_iterator(raw_mg.time);
            double error = abs(solution.position - ts->position);

            logger.logv("GDOP", solution.GDOP);
            logger.logv("Solution norm", abs(solution.position));
            logger.logv("True norm", abs(ts->position));
            logger.logv("Error", error);
            if (error > 10) logger.log("HIGH ERROR");
        }

        logger.log("");

        solution_states.push_back(solution);
    }
   
    error_type = '0';

    logger.log("Finished solving");
}

bool SatNav::check_raw(const RawMeasurement& raw_m) {
    if (!raw_m.is_present) {
        return false;
    }

    if (raw_m.qualflg != 0) {
        return false;
    }

    if (error_type != 'S') {
        double min_CN0 = std::min(raw_m.L1_CN0, raw_m.L2_CN0);
        double max_CN0 = std::max(raw_m.L1_CN0, raw_m.L2_CN0);

        if (min_CN0 < CN0_min_threshold || max_CN0 > CN0_max_threshold) {
            logger.log("Bad CN0: " + std::to_string(raw_m.prn_id) + ", " + std::to_string(min_CN0) + ", " + std::to_string(max_CN0));
            return false;
        }
    }

    if (error_type != 'F') {
        if (check_fading(raw_m)) {
            logger.logv("Fading", raw_m.prn_id);
            return false;
        }
    }

    return true;
}

RefinedMeasurement SatNav::refine_raw(const RawMeasurement& raw_m) {

    double pseudorange = raw_m.L1_range * C1 + raw_m.L2_range * C2;
    double carrier_phase = raw_m.L1_phase * C1 + raw_m.L2_phase * C2;

    double clock_error = handler.get_clock_error(raw_m.prn_id, raw_m.time);
    double relativistic_error = handler.get_relativistic_error(raw_m.prn_id, raw_m.time);

    pseudorange += clock_error * c + relativistic_error * c;
    double dt = pseudorange / c;

    State gps_state = handler.get_state(raw_m.prn_id, raw_m.time - dt);

    std::vector<double> gps_position = gps_state.position;
    std::vector<double> gps_velocity = gps_state.velocity;

    if (error_type != 's') {
        double phi = earth_rotation_rate * dt;
        Matrix rot = rotation(-phi, 'z');
        gps_position = rot * gps_position;
    }

    return {true, raw_m.time, raw_m.prn_id, pseudorange, carrier_phase, gps_position, gps_velocity};
    
}

    /*  Эффект                     t       t * v    t * c   
    (1) Задержка распространения - 80 мс - 240 км -  ---
    (2) Ошибка часов НКА         - 3 мс  - 10 м   - 1000 км
    (3) Релятивизм               - 20 нс - 60 мкм - 6 м

    При расчете ошибки часов НКА учитываем только (1)
    При расчете эфемерид учитываем (1), (2)
    При корректировке псевдодальности учитываем (1), (2), (3) */

SolutionState SatNav::calculate_solution(const RefinedMeasurementGroupped& ref_mg) const {
    SolutionState solution;
    solution.time = ref_mg.time;

    std::vector<double> PR;
    std::vector<std::vector<double>> X;
    size_t n = 0;

    for (const auto& ref_m : ref_mg.refined_measurements) {
        if (ref_m.is_present) {
            PR.push_back(ref_m.pseudorange);
            X.push_back(ref_m.gps_position);
            ++n;
        }
    }

    double eps = 0.5;
    std::vector<double> x = {0.0, 0.0, 0.0};
    double c_tau = 0.0;

    std::vector<double> U(n);
    Matrix B(n, 4);

    while (true) {
        for (size_t i = 0; i < n; ++i) {
            std::vector<double> DX = X[i] - x;
            double D = abs(DX);

            DX = DX / D;
            DX.push_back(1.0);

            U[i] = PR[i] - D;
            B.at(i) = DX;
        }

        Matrix B1(4, 4);
        try {
            B1 = (B.transpose() * B).inverse();
        } catch (const std::runtime_error&) {
            solution.is_solved = false;
            solution.failure_type = 'g';
            return solution;
        }

        double GDOP = sqrt(B1.trace());
        solution.GDOP = GDOP;

        if (GDOP > GDOP0) {
            solution.is_solved = false;
            solution.failure_type = 'G';
            return solution;
        }

        std::vector<double> dX = -B1 * B.transpose() * U;

        std::vector<double> dX_x = std::vector<double>(dX.begin(), dX.begin() + 3);
        double dX_c_tau = dX[3];

        double delta = abs(dX_x) + std::abs(dX_c_tau - c_tau);       
        if (delta < eps) {
            break;
        }

        x = x + dX_x;
        c_tau = dX_c_tau;
    }

    solution.position = x;
    solution.is_solved = true;

    return solution;
}

bool SatNav::check_fading(const RawMeasurement& raw_m) {
    unsigned t0 = raw_m.time;
    unsigned prn_index = raw_m.prn_id - 1;
    
    auto raw_mg_it = get_raw_measurement_groupped_iterator(t0);

    auto it_bwd = raw_mg_it;
    while(true) {
        if (it_bwd == raw_measurements_groupped.begin()) return true;
        it_bwd--;

        unsigned t = it_bwd->time;
        if (t0 - t > fadeout_threshold) break;
    
        if (!it_bwd->raw_measurements[prn_index].is_present) return true;
    }
    
    return false;
}

std::vector<unsigned> SatNav::check_low(const SolutionState& solution, const RefinedMeasurementGroupped& ref_mg) {
    std::vector<unsigned> low;

    for (unsigned prn_id = 1; prn_id <= 32; ++prn_id) {
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

const std::vector<State>& SatNav::get_true_states() const {
    return true_states;
}

const std::vector<SolutionState>& SatNav::get_solution_states() const {
    return solution_states;
}

const std::vector<RawMeasurementGroupped>& SatNav::get_raw_measurements_groupped() const {
    return raw_measurements_groupped;
}

const std::vector<RefinedMeasurementGroupped>& SatNav::get_refined_measurements_groupped() const {
    return refined_measurements_groupped;
}

const std::vector<AccelerationMeasurement>& SatNav::get_acceleration_measurements() const {
    return acceleration_measurements;
}

// сделать темплейт бинарного поиска
std::vector<State>::const_iterator SatNav::get_true_state_iterator(unsigned time) const {
    size_t n = true_states.size();

    std::ptrdiff_t low = 0;
    std::ptrdiff_t high = static_cast<std::ptrdiff_t>(n) - 1;
    while (low <= high) {
        std::ptrdiff_t mid = low + (high - low) / 2;

        if (true_states[static_cast<size_t>(mid)].time == time) {
            return true_states.cbegin() + mid;
        }

        if (true_states[static_cast<size_t>(mid)].time < time) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }

    return true_states.cend();
}

std::vector<SolutionState>::const_iterator SatNav::get_solution_state_iterator(unsigned time) const {
    size_t n = solution_states.size();

    std::ptrdiff_t low = 0;
    std::ptrdiff_t high = static_cast<std::ptrdiff_t>(n) - 1;
    while (low <= high) {
        std::ptrdiff_t mid = low + (high - low) / 2;

        if (solution_states[static_cast<size_t>(mid)].time == time) {
            return solution_states.cbegin() + mid;
        }

        if (solution_states[static_cast<size_t>(mid)].time < time) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }

    return solution_states.cend();
}

std::vector<RawMeasurementGroupped>::const_iterator SatNav::get_raw_measurement_groupped_iterator(unsigned time) const {
    size_t n = raw_measurements_groupped.size();

    std::ptrdiff_t low = 0;
    std::ptrdiff_t high = static_cast<std::ptrdiff_t>(n) - 1;
    while (low <= high) {
        std::ptrdiff_t mid = low + (high - low) / 2;

        if (raw_measurements_groupped[static_cast<size_t>(mid)].time == time) {
            return raw_measurements_groupped.cbegin() + mid;
        }

        if (raw_measurements_groupped[static_cast<size_t>(mid)].time < time) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }

    return raw_measurements_groupped.cend();
}
