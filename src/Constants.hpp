#pragma once

#include <Eigen/Dense>

namespace Constants {
    inline constexpr double c = 2.99792458e8;

    inline constexpr double mu = 3.986004188e14;
    inline constexpr double Re = 6378136.6;
    inline constexpr double omega = 7.292115e-5;
    inline constexpr double J2 = 1.08263e-3;

    inline constexpr double C1 = 2.54572778016;
    inline constexpr double C2 = -1.54572778016;
    inline constexpr double CN0_constant = 3.01029995664;

    inline const Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
    inline const Eigen::Matrix3d Omega = [] {
        Eigen::Matrix3d M;
        M << 0.0, -omega, 0.0, omega, 0.0, 0.0, 0.0, 0.0, 0.0;
        return M;
    }();

    inline const Eigen::Matrix3d Omega2 = [] {
        Eigen::Matrix3d M;
        double omega2 = -omega * omega;
        M << omega2, 0.0, 0.0, 0.0, omega2, 0.0, 0.0, 0.0, 0.0;
        return M;
    }();
}