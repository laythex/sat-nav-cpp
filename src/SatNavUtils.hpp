#pragma once

#include <Eigen/Dense>
#include <vector>

#include "Constants.hpp"

namespace SatNavUtils {
    bool check_time(unsigned t, unsigned ti, unsigned tf, unsigned t0);

    double ionofree(double l1, double l2);

    double calculate_delta(const Eigen::Vector3d& av, const Eigen::Vector3d& xv);

    Eigen::Matrix3d calculate_G(const Eigen::Vector3d& rv);

    Eigen::Matrix3d calculate_B();

    Eigen::VectorXd shifted_difference(const Eigen::VectorXd& v);
    Eigen::MatrixX3d shifted_difference(const Eigen::MatrixX3d& m);
}
