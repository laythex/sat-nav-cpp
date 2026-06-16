#pragma once

#include <Eigen/Dense>
#include <vector>

#include "SatNavUtils.hpp"
#include "Structures.hpp"

struct KalmanState {
    unsigned time = 0;

    Eigen::Vector3d R = Eigen::Vector3d::Zero();
    Eigen::Vector3d dR = Eigen::Vector3d::Zero();

    Eigen::Vector3d r_pas = Eigen::Vector3d::Zero();
    Eigen::Vector3d dr_pas = Eigen::Vector3d::Zero();

    Eigen::Vector3d r_act = Eigen::Vector3d::Zero();

    std::vector<bool> gps_mask = std::vector<bool>(32, false);
    std::vector<double> zeta1i_groupped = std::vector<double>(32, 0.0);
    std::vector<double> zeta2i_groupped = std::vector<double>(32, 0.0);
    std::vector<Eigen::Vector3d> Ci_groupped = std::vector<Eigen::Vector3d>(32, Eigen::Vector3d::Zero());

    Eigen::MatrixXd W = Eigen::MatrixXd::Zero(6, 6);

    void fill_measurements(const StandaloneData& data);

    DiscreteState get_dstate_rel() const;
    DiscreteState get_dstate_pas() const;
};