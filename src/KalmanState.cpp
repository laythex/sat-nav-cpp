#include "KalmanState.hpp"

void KalmanState::fill_measurements(const StandaloneData& data) {
    std::vector<size_t> present_prns;
    for (size_t prn_id = 1; prn_id <= 32; ++prn_id) {
        if (!data.is_fully_present_at(prn_id)) continue;
        present_prns.push_back(prn_id);
    }
    size_t n = present_prns.size();

    for (size_t i = 0; i < n; ++i) {
        size_t prn_index = present_prns[i] - 1;
        gps_mask[prn_index] = true;

        Eigen::Vector3d Sv = data.pas.position - data.pas.gps_states[prn_index].position;
        Eigen::Vector3d Xv = data.act.position - data.pas.position;
        Eigen::Vector3d Yv = data.act.gps_states[prn_index].position - data.pas.gps_states[prn_index].position;

        double delta = SatNavUtils::delta(Sv, Xv - Yv);

        double pr_difference = data.act.prs_L1[prn_index] - data.pas.prs_L1[prn_index];
        double cp_difference = data.act.cps_L1[prn_index] - data.pas.cps_L1[prn_index];

        Eigen::Vector3d Svn = Sv.normalized();

        zeta1i_groupped[prn_index] = pr_difference - delta + Svn.dot(Yv);
        zeta2i_groupped[prn_index] = cp_difference - delta + Svn.dot(Yv);
        Ci_groupped[prn_index] = Svn;
    }
}

DiscreteState KalmanState::get_dstate_rel() const { return {time, R, dR}; }

DiscreteState KalmanState::get_dstate_pas() const { return {time, r_pas, dr_pas}; }
