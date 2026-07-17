#include "SatNav.hpp"
#include "SatNavRel.hpp"

// Собирать полный путь либо здесь, либо в DataParser
SatNav::SatNav(const Date& date, const GRACE sat, const GPSHandler& handler)
    : date(date), handler(handler), logger("GPS") {
    bool is_fo = sat == GRACE::C || sat == GRACE::D;
    if (is_fo) {
        true_states = DataParser::load_grace_fo_gnv_data(date, sat);
        raw_measurements_groupped = DataParser::load_grace_fo_gps_data(date, sat);
    } else {
        true_states = DataParser::load_grace_gnv_data(date, sat);
        raw_measurements_groupped = DataParser::load_grace_gps_data(date, sat);
    }
}

SatNav::SatNav(const Date& date, const SWARM sat, const GPSHandler& handler)
    : date(date), handler(handler), logger("GPS") {
    true_states = DataParser::load_swarm_nav_data(date, sat);
    raw_measurements_groupped = DataParser::load_swarm_gps_data(date, sat);
}

SatNav::SatNav(const SatNav& sn)
    : date(sn.date), handler(sn.handler), logger("GPS"), true_states(sn.true_states),
      raw_measurements_groupped(sn.raw_measurements_groupped), acceleration_measurements(sn.acceleration_measurements) {
}

void SatNav::solve(IntendedError error) {
    this->error = error;

    logger.logv("Starting GPS", date);
    logger.lnbr();

    unsigned t0 = raw_measurements_groupped.front().time;
    for (const RawMeasurementGroupped& raw_mg : raw_measurements_groupped) {
        if (raw_mg.time - t0 < Config::ti || raw_mg.time - t0 >= Config::tf) continue;

        logger.logv("Time", raw_mg.time);

        SolutionState solution = find_solution(raw_mg);
        solution_states.push_back(solution);

        State ts = *get_true_state_iterator(raw_mg.time);
        logger.logv("Status", solution.status);
        logger.logv("GDOP", solution.GDOP);
        logger.logv("Error", (solution.position - ts.position).norm());
        logger.lnbr();
    }

    logger.logv("Finishing GPS", date);
}

SolutionState SatNav::find_solution(const RawMeasurementGroupped& raw_mg) {
    std::vector<size_t> present_prns = check_raw_groupped(raw_mg);
    RefinedMeasurementGroupped ref_mg = refine_raw_groupped(raw_mg, present_prns);
    refined_measurements_groupped.push_back(ref_mg);

    size_t n = present_prns.size();
    logger.logv("PRN count", n);

    double GDOP;
    double dt_ASN = 0.0;
    Eigen::Vector3d X = Eigen::Vector3d::Zero();
    double residual = -1.0;

    Eigen::Index ne = static_cast<Eigen::Index>(n);
    Eigen::VectorXd U(ne);
    Eigen::MatrixX4d B(ne, 4);

    while (true) {
        for (size_t i = 0; i < n; ++i) {
            size_t prn_index = present_prns[i] - 1;

            const RefinedMeasurement& ref_m = ref_mg.measurements[prn_index];

            Eigen::Vector3d Xi = ref_m.gps_state.position;
            if (error != IntendedError::SAGNAC) {
                Xi -= Constants::Omega * (ref_m.toa_ASN - ref_m.tot - dt_ASN) * ref_m.gps_state.position;
            }

            Eigen::Index ie = static_cast<Eigen::Index>(i);
            U(ie) = ref_m.pr_refined - (Xi - X).norm();
            B.row(ie) = (Xi - X).normalized().homogeneous();
        }

        Eigen::FullPivLU<Eigen::Matrix4d> lu(B.transpose() * B);
        if (!lu.isInvertible()) {
            return {{raw_mg.time}, SolutionStatus::SINGULAR_B1};
        }

        Eigen::Matrix4d B1 = lu.inverse();

        GDOP = sqrt(B1.trace());
        if (error != IntendedError::GDOP) {
            if (GDOP > Config::GDOP0) {
                return {{raw_mg.time}, SolutionStatus::HIGH_GDOP, GDOP};
            }
        }

        Eigen::Vector4d dX = -B1 * B.transpose() * U;
        Eigen::Vector3d dX_X = dX.head<3>();
        double dX_dt_ASN = dX[3] / Constants::c;

        X += dX_X;
        dt_ASN = dX_dt_ASN;

        Eigen::MatrixXd In = Eigen::MatrixXd::Identity(ne, ne);
        residual = ((In - B * B1 * B.transpose()) * U).norm();

        double delta = dX_X.norm() + std::abs(dX_dt_ASN - dt_ASN) * Constants::c;
        if (delta < Config::solution_tolerance) break;
    }

    if (residual > Config::residual_threshold && error == IntendedError::NONE) {
        return {{raw_mg.time}, SolutionStatus::DEFAULT};
    }

    SolutionState ss = {{raw_mg.time, X}, SolutionStatus::OK, GDOP, dt_ASN * Constants::c, residual};

    for (size_t i = 0; i < n; ++i) {
        size_t prn_index = present_prns[i] - 1;

        const RefinedMeasurement& ref_m = ref_mg.measurements[prn_index];

        Eigen::Vector3d Xi = ref_m.gps_state.position;
        if (error != IntendedError::SAGNAC) {
            Xi -= Constants::Omega * (ref_m.toa_ASN - ref_m.tot - dt_ASN) * ref_m.gps_state.position;
        }

        ss.gps_mask[prn_index] = true;
        ss.gps_states[prn_index] = {raw_mg.time, Xi};
        ss.prs_L1[prn_index] = raw_mg.measurements[prn_index].L1_range;
        ss.cps_L1[prn_index] = raw_mg.measurements[prn_index].L1_phase;
    }

    return ss;
}

std::vector<size_t> SatNav::check_raw_groupped(const RawMeasurementGroupped& raw_mg) {
    std::vector<size_t> present_prns;

    for (size_t prn_id = 1; prn_id <= 32; ++prn_id) {
        size_t prn_index = prn_id - 1;
        const RawMeasurement& raw_m = raw_mg.measurements[prn_index];

        MeasurementStatus m_status = check_raw(raw_m);
        if (m_status != MeasurementStatus::NOT_PRESENT) {
            logger.logv("Checking " + std::to_string(prn_id), m_status);
        }
        if (m_status != MeasurementStatus::OK) continue;

        present_prns.push_back(prn_id);
    }

    return present_prns;
}

MeasurementStatus SatNav::check_raw(const RawMeasurement& raw_m) const {
    if (!raw_m.is_present) {
        return MeasurementStatus::NOT_PRESENT;
    }

    if (error != IntendedError::QUALITY_FLAG) {
        if (raw_m.qualflg != 0) {
            return MeasurementStatus::QUALITY_FLAG;
        }
    }

    if (error != IntendedError::NOISY) {
        double CN0_min = std::min(raw_m.L1_CN0, raw_m.L2_CN0);
        double CN0_max = std::max(raw_m.L1_CN0, raw_m.L2_CN0);
        if (CN0_min < Config::CN0_min_threshold || CN0_max > Config::CN0_max_threshold) {
            return MeasurementStatus::CN0;
        }

        if (check_fadein(raw_m)) {
            return MeasurementStatus::FADE_IN;
        }
    }

    return MeasurementStatus::OK;
}

RefinedMeasurementGroupped SatNav::refine_raw_groupped(const RawMeasurementGroupped& raw_mg,
                                                       const std::vector<size_t>& present_prns) const {
    RefinedMeasurementGroupped ref_mg = {raw_mg.time};

    size_t n = present_prns.size();
    for (size_t i = 0; i < n; ++i) {
        size_t prn_index = present_prns[i] - 1;

        RefinedMeasurement ref_m = refine_raw(raw_mg.measurements[prn_index]);

        if (error != IntendedError::HATCH) {
            ref_m = hatch_filter(ref_m);
        }

        ref_mg.measurements[prn_index] = ref_m;
    }

    return ref_mg;
}

RefinedMeasurement SatNav::refine_raw(const RawMeasurement& raw_m) const {
    unsigned prn_id = raw_m.prn_id;
    unsigned toa_ASN = raw_m.time;

    double clock_error = handler.get_clock_error(prn_id, toa_ASN);
    double relativistic_error = handler.get_relativistic_error(prn_id, toa_ASN);

    if (error == IntendedError::GPS_CLOCK) {
        clock_error = 0.0;
    }

    if (error == IntendedError::RELATIVISTIC) {
        relativistic_error = 0.0;
    }

    double tot = toa_ASN - raw_m.L1_range / Constants::c - clock_error;

    double pr_ionofree = SatNavUtils::ionofree(raw_m.L1_range, raw_m.L2_range);
    double cp_ionofree = SatNavUtils::ionofree(raw_m.L1_phase, raw_m.L2_phase);

    if (error == IntendedError::IONOSPHERIC) {
        pr_ionofree = raw_m.L1_range;
        cp_ionofree = raw_m.L1_phase;
    }

    double pr_refined = pr_ionofree + (clock_error - relativistic_error) * Constants::c;

    State gps_state = handler.get_state(prn_id, tot);

    return {MeasurementStatus::OK, prn_id, toa_ASN, tot, pr_refined, cp_ionofree, gps_state};
}

bool SatNav::check_fadein(const RawMeasurement& raw_m) const {
    unsigned t0 = raw_m.time;
    unsigned prn_index = raw_m.prn_id - 1;

    auto raw_mg_it = get_raw_measurement_groupped_iterator(t0);

    auto it_bwd = raw_mg_it;
    while (true) {
        if (it_bwd == raw_measurements_groupped.begin()) return true;
        it_bwd--;

        unsigned t = it_bwd->time;
        if (t0 - t > Config::fadein_threshold) break;

        if (!it_bwd->measurements[prn_index].is_present) return true;
    }

    return false;
}

RefinedMeasurement SatNav::hatch_filter(const RefinedMeasurement& ref_m) const {
    unsigned prn_index = ref_m.prn_id - 1;

    if (refined_measurements_groupped.empty()) return ref_m;
    const RefinedMeasurement& ref_m_prev = refined_measurements_groupped.back().measurements[prn_index];
    if (ref_m_prev.status != MeasurementStatus::OK) return ref_m;

    double pr_prev = ref_m_prev.pr_refined;
    double cp_prev = ref_m_prev.cp_refined;
    double cp_delta = ref_m.cp_refined - cp_prev;
    double pr_hatch = Config::hatch_constant * ref_m.pr_refined + (1.0 - Config::hatch_constant) * (pr_prev + cp_delta);

    RefinedMeasurement ref_m_hatch = ref_m;
    ref_m_hatch.pr_refined = pr_hatch;
    return ref_m_hatch;
}

const std::vector<State>& SatNav::get_true_states() const { return true_states; }

const std::vector<SolutionState>& SatNav::get_solution_states() const { return solution_states; }

const std::vector<RawMeasurementGroupped>& SatNav::get_raw_measurements_groupped() const {
    return raw_measurements_groupped;
}

const std::vector<RefinedMeasurementGroupped>& SatNav::get_refined_measurements_groupped() const {
    return refined_measurements_groupped;
}

const std::vector<AccelerationMeasurement>& SatNav::get_acceleration_measurements() const {
    return acceleration_measurements;
}

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

Date SatNav::get_date() const { return date; }
