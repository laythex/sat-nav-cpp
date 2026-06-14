#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "Constants.hpp"
#include "Structures.hpp"

namespace SatNavUtils {
    double ionofree(double l1, double l2);

    double delta(const Eigen::Vector3d& av, const Eigen::Vector3d& xv);

    Eigen::Vector3d acceleration(const Eigen::Vector3d& rv);

    DiscreteState step_inc(const DiscreteState& s, double dt);
    DiscreteState step_inc_rel(const DiscreteState& s_rel, const Eigen::Vector3d rv_pas, double dt);

    Eigen::Matrix3d G(const Eigen::Vector3d& rv);
    Eigen::MatrixXd iF(const Eigen::Matrix3d& G, double dt);

    Eigen::VectorXd shifted_difference(const Eigen::VectorXd& v);
    Eigen::MatrixX3d shifted_difference(const Eigen::MatrixX3d& m);

    Eigen::VectorXd form_groupped(const std::vector<double>& groupped, const std::vector<size_t>& present_prns);
    Eigen::MatrixX3d form_groupped(const std::vector<Eigen::Vector3d>& groupped,
                                   const std::vector<size_t>& present_prns);
}
