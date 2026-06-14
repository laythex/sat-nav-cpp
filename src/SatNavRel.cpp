#include "SatNavRel.hpp"
#include "SatNav.hpp"

SatNavRel::SatNavRel(SatNav& passive, SatNav& active) : pas(passive), act(active), logger("RGPS") {}

SatNavRel::SatNavRel(const SatNavRel& sn) : pas(sn.pas), act(sn.act), true_states(sn.true_states), logger("RGPS") {}

void SatNavRel::solve_standalones(unsigned ti, unsigned tf) {
    pas.solve(ti, tf);
    act.solve(ti, tf);
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

void SatNavRel::solve(unsigned ti, unsigned tf) {
    solve_standalones(ti, tf);
    form_true_states();
    form_standalone_data();

    logger.logv("Starting RGPS", pas.date);
    logger.lnbr();

    unsigned t0 = standalone_datas[0].time;
    for (const StandaloneData& data : standalone_datas) {
        if (!SatNavUtils::check_time(data.time, ti, tf, t0)) continue;

        logger.logv("Time", data.time);

        SolutionStateRel solution = find_solution_ls(data);
        solution_states.push_back(solution);

        State ts = *get_true_state_iterator(data.time);
        logger.logv("Status", solution.status);
        logger.logv("Error", (solution.position - ts.position).norm());
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

        double delta = SatNavUtils::calculate_delta(Sv, Xv - Yv);

        double pr_pas = data.pas.prs_L1[prn_index];
        double pr_act = data.act.prs_L1[prn_index];
        double pr_difference = pr_act - pr_pas;

        Eigen::Vector3d Svn = Sv.normalized();

        Eigen::Index ei = static_cast<Eigen::Index>(i);
        Ui[ei] = pr_difference - delta + Svn.dot(Yv);
        Ci.row(ei) = Svn;
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