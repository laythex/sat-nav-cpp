#include "SatNav.hpp"
#include "SatNavRel.hpp"

// Собирать полный путь либо здесь, либо в DataParser
SatNav::SatNav(const Date& date, const GRACE sat, const GPSHandler& handler) : handler(handler), logger("GPS") {
    bool is_fo = sat == GRACE::C || sat == GRACE::D;
    if (is_fo) {
        true_states = DataParser::load_grace_fo_gnv_data(date, sat);
        raw_measurements_groupped = DataParser::load_grace_fo_gps_data(date, sat);
    } else {
        true_states = DataParser::load_grace_gnv_data(date, sat);
        raw_measurements_groupped = DataParser::load_grace_gps_data(date, sat);
    }
}

SatNav::SatNav(const Date& date, const SWARM sat, const GPSHandler& handler) : handler(handler), logger("GPS") {
    true_states = DataParser::load_swarm_nav_data(date, sat);
    raw_measurements_groupped = DataParser::load_swarm_gps_data(date, sat);
}

SatNav::SatNav(const SatNav& sn) : handler(sn.handler), 
                                   true_states(sn.true_states), 
                                   raw_measurements_groupped(sn.raw_measurements_groupped),
                                   acceleration_measurements(sn.acceleration_measurements), logger("GPS") { }

void SatNav::solve(unsigned ti, unsigned tf, IntendedError error) {
    logger.log("Starting GPS solving");

    this->error = error;

    unsigned t0 = raw_measurements_groupped[0].time;
    for (const auto& raw_mg : raw_measurements_groupped) {
        if (ti > 0 || tf > 0) {
            if ((raw_mg.time - t0) < ti || (raw_mg.time - t0) >= tf) continue;
        }

        logger.logv("Current time", raw_mg.time);

        std::vector<RefinedMeasurement> ref_ms(32);
        for (unsigned prn_id = 1; prn_id <= 32; ++prn_id) {
            unsigned prn_index = prn_id - 1;

            RefinedMeasurement ref_m = refine_raw(raw_mg.raw_measurements[prn_index]);

            if (ref_m.status != MeasurementStatus::NOT_PRESENT) {
                logger.log("Refining at " + std::to_string(prn_id) + ": " + to_string(ref_m.status));
            }

            if (ref_m.status == MeasurementStatus::OK) {
                ref_m = hatch_filter(ref_m);
            }

            ref_ms[prn_index] = ref_m;
        }

        RefinedMeasurementGroupped ref_mg = {raw_mg.time, ref_ms};
        refined_measurements_groupped.push_back(ref_mg);

        SolutionState solution = calculate_solution(ref_mg);

        // if (solution.status == SolutionStatus::OK) {
        //     std::vector<unsigned> low = check_low(solution, ref_mg);
        //     if (low.size() > 0) {
        //         for (unsigned prn_id : low) {
        //             unsigned prn_index = prn_id - 1;
        //             ref_mg.refined_measurements[prn_index].status = MeasurementStatus::FADING;
        //         }
        //         solution = calculate_solution(ref_mg);
        //     }
        // }

        logger.logv("Solution status", solution.status);

        if (solution.status == SolutionStatus::OK) {
            auto ts = get_true_state_iterator(raw_mg.time);
            double error = abs(solution.position - ts->position);

            logger.logv("GDOP", solution.GDOP);
            logger.logv("Solution norm", abs(solution.position));
            logger.logv("True norm", abs(ts->position));
            logger.logv("Error", error);
            if (error > 10) logger.log("HIGH ERROR");
        }

        logger.lnbr();

        solution_states.push_back(solution);
    }
   
    // error_type = '0';

    logger.log("Finished GPS solving");
}

MeasurementStatus SatNav::check_raw(const RawMeasurement& raw_m) {
    if (!raw_m.is_present) {
        return MeasurementStatus::NOT_PRESENT;
    }

    if (raw_m.qualflg != 0) {
        return MeasurementStatus::QUALITY_FLAG;
    }

    double min_CN0 = std::min(raw_m.L1_CN0, raw_m.L2_CN0);
    double max_CN0 = std::max(raw_m.L1_CN0, raw_m.L2_CN0);
    if (min_CN0 < CN0_min_threshold || max_CN0 > CN0_max_threshold) {
        return MeasurementStatus::CN0;
    }

    // if (check_fading(raw_m)) {
    //     logger.logv("Fading at ", raw_m.prn_id);
    //     return MeasurementStatus::FADING;
    // }

    return MeasurementStatus::OK;
}

RefinedMeasurement SatNav::refine_raw(const RawMeasurement& raw_m) {
    MeasurementStatus status = check_raw(raw_m);
    if (status != MeasurementStatus::OK) {
        return {.status = status};
    }

    unsigned toa_ASN = raw_m.time;

    double clock_error = handler.get_clock_error(raw_m.prn_id, toa_ASN);
    double relativistic_error = handler.get_relativistic_error(raw_m.prn_id, toa_ASN);

    double tot = toa_ASN - raw_m.L1_range / c - clock_error;

    double pr_ionofree = raw_m.L1_range * C1 + raw_m.L2_range * C2;
    double cp_ionofree = raw_m.L1_phase * C1 + raw_m.L2_phase * C2;

    double pr_refined = pr_ionofree + (clock_error - relativistic_error) * c;

    State gps_state = handler.get_state(raw_m.prn_id, tot);
    std::vector<double> gps_position = gps_state.position;
    std::vector<double> gps_velocity = gps_state.velocity;
    
    return {status, raw_m.prn_id, toa_ASN, tot, pr_refined, cp_ionofree, gps_position, gps_velocity};    
}

SolutionState SatNav::calculate_solution(RefinedMeasurementGroupped& ref_mg) const {
    SolutionState solution;
    solution.time = ref_mg.time;

    double dt_ASN = 0.0;
    std::vector<double> X(3, 0.0);

    while (true) {
        std::vector<double> U(0);
        std::vector<std::vector<double>> B_v;

        for (const auto& ref_m : ref_mg.refined_measurements) {
            if (ref_m.status != MeasurementStatus::OK) continue;

            std::vector<double> X_i = (E3 + Omega * (ref_m.toa_ASN - ref_m.tot - dt_ASN)) * ref_m.gps_position;

            std::vector<double> a = normalize(X_i - X);
            double b = abs(X_i - X);

            U.push_back(ref_m.pseudorange - b);
            B_v.push_back({a[0], a[1], a[2], 1});
        }

        Matrix B(B_v);
        Matrix B1(4, 4);
        try {
            B1 = (B.transpose() * B).inverse();
        } catch (const std::runtime_error&) {
            solution.status = SolutionStatus::SINGULAR_B1;
            return solution;
        }

        double GDOP = sqrt(B1.trace());
        solution.GDOP = GDOP;

        if (GDOP > GDOP0 || GDOP == 0) {
            solution.status = SolutionStatus::HIGH_GDOP;
            return solution;
        }

        std::vector<double> dX = -B1 * B.transpose() * U;
        std::vector<double> dX_X = {dX[0], dX[1], dX[2]};
        double dX_dt_ASN = dX[3] / c;

        double delta = abs(dX_X) + abs(dX_dt_ASN - dt_ASN) * c;
        if (delta < solution_tolerance) {
            break;
        }

        X = X + dX_X;
        dt_ASN = dX_dt_ASN;
    }

    solution.position = X;
    solution.status = SolutionStatus::OK;

    return solution;
}

// bool SatNav::check_fading(const RawMeasurement& raw_m) {
//     unsigned t0 = raw_m.time;
//     unsigned prn_index = raw_m.prn_id - 1;
    
//     auto raw_mg_it = get_raw_measurement_groupped_iterator(t0);

//     auto it_bwd = raw_mg_it;
//     while(true) {
//         if (it_bwd == raw_measurements_groupped.begin()) return true;
//         it_bwd--;

//         unsigned t = it_bwd->time;
//         if (t0 - t > fadeout_threshold) break;
    
//         if (!it_bwd->raw_measurements[prn_index].is_present) return true;
//     }
    
//     return false;
// }

// std::vector<unsigned> SatNav::check_low(const SolutionState& solution, const RefinedMeasurementGroupped& ref_mg) {
//     std::vector<unsigned> low;

//     for (unsigned prn_id = 1; prn_id <= 32; ++prn_id) {
//         unsigned prn_index = prn_id - 1;
//         RefinedMeasurement ref_m = ref_mg.refined_measurements[prn_index];

//         if (ref_m.status != ) continue;
        
//         std::vector<double> gps_relative = ref_m.gps_position - solution.position;
//         double zenith_angle = angle(solution.position, gps_relative) * 180.0 * M_1_PI;
//         if (90.0 - zenith_angle < mask_angle) {
//             low.push_back(prn_id);
//         }
//     }

//     return low;
// }

RefinedMeasurement SatNav::hatch_filter(const RefinedMeasurement& ref_m) {
    unsigned prn_index = ref_m.prn_id - 1;

    RefinedMeasurement ref_m_filtered = ref_m;

    size_t n = refined_measurements_groupped.size();
    if (n < 2) {
        return ref_m_filtered;
    }

    RefinedMeasurement ref_m_previous = refined_measurements_groupped[n - 2].refined_measurements[prn_index];  
    if (ref_m_previous.status != MeasurementStatus::OK) {
        return ref_m_filtered;
    }

    double pseudorange_previous = ref_m_previous.pseudorange;
    double carrier_phase_previous = ref_m_previous.carrier_phase;
    double pseudorande_filtered = pseudorange_previous + ref_m.carrier_phase - carrier_phase_previous;
    ref_m_filtered.pseudorange = (1.0 / hatch_constant) * ref_m.pseudorange + 
                                 (1.0 - 1.0 / hatch_constant) * pseudorande_filtered;
    
    return ref_m_filtered;
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
