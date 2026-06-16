#include "SatNavRel.hpp"
#include "SatNav.hpp"

SatNavRel::SatNavRel(SatNav& passive, SatNav& active) : pas(passive), act(active), logger("RGPS") {}

SatNavRel::SatNavRel(const SatNavRel& sn) : pas(sn.pas), act(sn.act), true_states(sn.true_states), logger("RGPS") {}

void SatNavRel::solve_standalones() {
    pas.solve();
    act.solve();
}

void SatNavRel::form_true_states() {
    for (const State& ts_pas : pas.get_true_states()) {
        unsigned time = ts_pas.time;

        auto ts_act_it = act.get_true_state_iterator(time);
        if (ts_act_it == act.get_true_states().end()) continue;
        State ts_act = *ts_act_it;

        true_states.push_back({time, ts_act.position - ts_pas.position, ts_act.velocity - ts_pas.velocity});
    }
}

void SatNavRel::form_standalone_data() {
    for (const SolutionState& ss_pas : pas.get_solution_states()) {
        unsigned time = ss_pas.time;

        auto ss_act_it = act.get_solution_state_iterator(time);
        if (ss_act_it == act.get_solution_states().end()) continue;
        SolutionState ss_act = *ss_act_it;

        standalone_datas.push_back({time, ss_pas, ss_act});
    }
}

void SatNavRel::solve() {
    solve_standalones();
    form_true_states();
    form_standalone_data();

    logger.logv("Starting RGPS", pas.date);
    logger.lnbr();

    for (const StandaloneData& data : standalone_datas) {
        logger.logv("Time", data.time);

        SolutionStateRel solution = Config::use_kf ? find_solution_kf(data) : find_solution_ls(data);
        solution_states.push_back(solution);

        State ts = *get_true_state_iterator(data.time);
        logger.logv("Status", solution.status);
        logger.logv("Error", (solution.position - ts.position).norm());
        if ((solution.position - ts.position).norm() > 1.5 && solution.status == SolutionStatusRel::OK) {
            logger.log("High error");
        }
        logger.lnbr();
    }

    logger.logv("Finishing RGPS", pas.date);
}

SolutionStateRel SatNavRel::find_solution_ls(const StandaloneData& data) {
    if (!data.is_fully_solved()) {
        return {{data.time}, SolutionStatusRel::NO_STANDALONE};
    }

    std::vector<size_t> present_prns;
    for (size_t prn_id = 1; prn_id <= 32; ++prn_id) {
        if (!data.is_fully_present_at(prn_id)) continue;
        present_prns.push_back(prn_id);
    }
    size_t n = present_prns.size();

    Eigen::VectorXd Ui(n);
    Eigen::MatrixX3d Ci(n, 3);

    for (size_t i = 0; i < n; ++i) {
        size_t prn_index = present_prns[i] - 1;

        Eigen::Vector3d Sv = data.pas.position - data.pas.gps_states[prn_index].position;
        Eigen::Vector3d Xv = data.act.position - data.pas.position;
        Eigen::Vector3d Yv = data.act.gps_states[prn_index].position - data.pas.gps_states[prn_index].position;

        double delta = SatNavUtils::delta(Sv, Xv - Yv);

        double pr_pas = data.pas.prs_L1[prn_index];
        double pr_act = data.act.prs_L1[prn_index];
        double pr_difference = pr_act - pr_pas;
        Eigen::Vector3d Svn = Sv.normalized();

        Eigen::Index ie = static_cast<Eigen::Index>(i);
        Ui(ie) = pr_difference - delta + Svn.dot(Yv);
        Ci.row(ie) = Svn;
    }

    Eigen::VectorXd U = SatNavUtils::shifted_difference(Ui);
    Eigen::MatrixX3d C = SatNavUtils::shifted_difference(Ci);

    Eigen::FullPivLU<Eigen::Matrix3d> lu(C.transpose() * C);
    if (!lu.isInvertible()) {
        return {{data.time}, SolutionStatusRel::SINGULAR_C1};
    }

    Eigen::Matrix3d C1 = lu.inverse();
    Eigen::Vector3d X = C1 * C.transpose() * U;

    return {{data.time, X}, SolutionStatusRel::OK};
}

SolutionStateRel SatNavRel::find_solution_kf(const StandaloneData& data) {
    if (kalman_states.size() < 2) {
        if (!data.is_fully_solved()) {
            logger.log("Cannot initialize");
            return {{data.time}, SolutionStatusRel::NO_STANDALONE};
        }

        logger.log("Initializing");
        kalman_states.push_back(initialize_kf(data));
        return {{data.time, data.act.position - data.pas.position}, SolutionStatusRel::INITIALIZING};
    }

    KalmanState ks;
    const KalmanState& ks_prev = kalman_states.back();
    ks.time = data.time;
    ks.fill_measurements(data);

    double res_max = std::max(data.pas.residual, data.act.residual);
    logger.logv("res max", res_max);
    Eigen::MatrixXd lambda(6, 6);
    lambda << Constants::I3 * Config::lambda_r, Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Zero(),
        Constants::I3 * Config::lambda_v;

    SolutionStatusRel status;
    if (!data.is_fully_solved() || res_max > Config::residual_threshold) {
        DiscreteState s_pas = SatNavUtils::step_inc(ks_prev.get_dstate_pas(), ks.time - ks_prev.time);
        ks.r_pas = s_pas.position;
        ks.dr_pas = s_pas.increment;
        status = SolutionStatusRel::MODEL_STEP;
        lambda *= Config::lambda_model;
    } else {
        ks.r_pas = data.pas.position;
        ks.dr_pas = ks.r_pas - ks_prev.r_pas;
        status = SolutionStatusRel::OK;
    }

    Eigen::VectorXd xi_hat = form_state(ks_prev, ks.time - ks_prev.time);
    Eigen::Matrix3d G = SatNavUtils::G(ks_prev.r_pas);
    Eigen::MatrixXd iF = SatNavUtils::iF(G, ks.time - ks_prev.time);

    auto measurements = form_measurements(ks, ks_prev);
    Eigen::VectorXd zeta = measurements.first;
    Eigen::MatrixXd D = measurements.second;
    Eigen::VectorXd zeta_hat = D * xi_hat;
    Eigen::Index m2 = zeta.size();

    Eigen::FullPivLU<Eigen::MatrixXd> lu(D.transpose() * D);
    if (!lu.isInvertible()) {
        status = SolutionStatusRel::MODEL_STEP;
    } else {
        Eigen::MatrixXd D1 = lu.inverse();
        Eigen::MatrixXd Im2 = Eigen::MatrixXd::Identity(m2, m2);
        double res_xi = ((Im2 - D * D1 * D.transpose()) * zeta).norm();
        double res_zeta = res_xi;
        logger.logv("res zeta", res_zeta);
        logger.logv("trace", D1.trace());
        if (res_zeta > Config::residual_rel_threshold || D1.trace() > Config::trace_threshold) {
            status = SolutionStatusRel::MODEL_STEP;
        }
    }

    double ir = 0.0;
    if (status == SolutionStatusRel::OK) {
        ir = 1.0;
    }

    ks.W = iF.transpose() * lambda * ks_prev.W * lambda * iF + D.transpose() * D * ir;
    Eigen::VectorXd xi_star = xi_hat + ks.W.inverse() * D.transpose() * (zeta - zeta_hat) * ir;

    ks.R = xi_star.head(3);
    ks.dR = xi_star.tail(3);

    logger.logv("m", m2 / 2);

    kalman_states.push_back(ks);

    return {{data.time, ks.R}, status};
}

KalmanState SatNavRel::initialize_kf(const StandaloneData& data) const {
    KalmanState ks;
    ks.time = data.time;
    ks.r_pas = data.pas.position;
    ks.R = data.act.position - data.pas.position;

    if (kalman_states.empty()) return ks;

    const KalmanState& ks_prev = kalman_states.back();
    ks.dr_pas = ks.r_pas - ks_prev.r_pas;
    ks.dR = ks.R - ks_prev.R;
    ks.fill_measurements(data);

    return ks;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> SatNavRel::form_measurements(const KalmanState& ks,
                                                                         const KalmanState& ks_prev) const {
    std::vector<size_t> common_prns;
    for (size_t prn_id = 1; prn_id <= 32; ++prn_id) {
        size_t prn_index = prn_id - 1;
        if (!ks.gps_mask[prn_index] || !ks_prev.gps_mask[prn_index]) continue;
        common_prns.push_back(prn_id);
    }
    size_t m = common_prns.size();

    Eigen::VectorXd zeta1 =
        SatNavUtils::shifted_difference(SatNavUtils::form_groupped(ks.zeta1i_groupped, common_prns));
    Eigen::VectorXd zeta2 =
        SatNavUtils::shifted_difference(SatNavUtils::form_groupped(ks.zeta2i_groupped, common_prns));
    Eigen::VectorXd zeta2_prev =
        SatNavUtils::shifted_difference(SatNavUtils::form_groupped(ks_prev.zeta2i_groupped, common_prns));
    Eigen::MatrixX3d C = SatNavUtils::shifted_difference(SatNavUtils::form_groupped(ks.Ci_groupped, common_prns));
    Eigen::MatrixX3d C_prev =
        SatNavUtils::shifted_difference(SatNavUtils::form_groupped(ks_prev.Ci_groupped, common_prns));

    Eigen::Index me = static_cast<Eigen::Index>(m);
    Eigen::VectorXd zeta(me * 2);
    Eigen::MatrixXd D(me * 2, 6);
    zeta << zeta1, zeta2 - zeta2_prev - (C - C_prev) * ks_prev.R;
    D << C, Eigen::MatrixX3d::Zero(me, 3), Eigen::MatrixX3d::Zero(me, 3), C;

    return {zeta, D};
}

Eigen::VectorXd SatNavRel::form_state(const KalmanState& ks_prev, double dt) const {
    State s_pas_ph = {ks_prev.time - 5, ks_prev.r_pas - ks_prev.dr_pas * 0.5, ks_prev.dr_pas / dt};
    State s_rel_ph = {ks_prev.time - 5, ks_prev.R - ks_prev.dR * 0.5, ks_prev.dR / dt};

    Eigen::Vector3d vv_prev = SatNavUtils::step_rk4(s_pas_ph, 5.0).velocity;
    Eigen::Vector3d Vv_prev = SatNavUtils::step_rk4_rel(s_rel_ph, s_pas_ph, 5.0).velocity;

    State s_pas_prev = {ks_prev.time, ks_prev.r_pas, vv_prev};
    State s_rel_prev = {ks_prev.time, ks_prev.R, Vv_prev};

    State s_rel_h = SatNavUtils::step_rk4_rel(s_rel_prev, s_pas_prev, dt * 0.5);
    State s_rel = SatNavUtils::step_rk4_rel(s_rel_prev, s_pas_prev, dt);

    Eigen::VectorXd xi(6);
    xi << s_rel.position, s_rel_h.velocity * dt;

    return xi;
}

const std::vector<State>& SatNavRel::get_true_states() const { return true_states; }

const std::vector<SolutionStateRel>& SatNavRel::get_solution_states() const { return solution_states; }

std::vector<State>::const_iterator SatNavRel::get_true_state_iterator(unsigned time) const {
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