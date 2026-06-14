#pragma once

#include <Eigen/Dense>
#include <vector>

#include "Constants.hpp"
#include "Structures.hpp"

class Propagator {
public:
    Propagator(const State& initial, bool useJ2 = true);

    State step_rk4(double step);
    State step_inc(double step);

private:
    Eigen::Vector3d acceleration(const Eigen::Vector3d& rv) const;
    Eigen::VectorXd rhs_rk4(const Eigen::VectorXd& sv) const;

    State initial;
    bool useJ2;
};