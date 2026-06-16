#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "DataParser.hpp"
#include "GPSHandler.hpp"
#include "KalmanState.hpp"
#include "Logger.hpp"
#include "SatNav.hpp"
#include "SatNavUtils.hpp"
#include "Structures.hpp"

class SatNav;

class SatNavRel {

public:
    SatNavRel(SatNav& passive, SatNav& active);
    SatNavRel(const SatNavRel& sn);

    void solve();

    const std::vector<State>& get_true_states() const;
    const std::vector<SolutionStateRel>& get_solution_states() const;

    std::vector<State>::const_iterator get_true_state_iterator(unsigned time) const;

    // Перенести в приват
    SatNav& pas;
    SatNav& act;

private:
    enum class SatType { PASSIVE, ACTIVE };

    void solve_standalones();
    void form_true_states();
    void form_standalone_data();

    SolutionStateRel find_solution_ls(const StandaloneData& data);
    SolutionStateRel find_solution_kf(const StandaloneData& data);

    KalmanState initialize_kf(const StandaloneData& data) const;

    std::pair<Eigen::VectorXd, Eigen::MatrixXd> form_measurements(const KalmanState& ks,
                                                                  const KalmanState& ks_prev) const;
    Eigen::VectorXd form_state(const KalmanState& ks_prev, double dt) const;

    std::vector<State> true_states;
    std::vector<StandaloneData> standalone_datas;
    std::vector<SolutionStateRel> solution_states;
    std::vector<KalmanState> kalman_states;

    Logger logger;
};
