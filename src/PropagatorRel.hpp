#pragma once

#include <Eigen/Dense>
#include <vector>

#include "Constants.hpp"
#include "Propagator.hpp"
#include "Structures.hpp"

class PropagatorRel {
public:
    PropagatorRel(const State& initial_rel);

    State step_rk4_rel(double step, const State& s_pas);
    State step_inc_rel(double step, const State& s_pas);

private:
    Eigen::Vector3d acceleration(const Eigen::Vector3d& rv) const;
    Eigen::Vector3d acceleration_rel(const Eigen::Vector3d& Rv, const Eigen::Vector3d& rv_pas) const;
    Eigen::VectorXd rhs_rk4_rel(const Eigen::VectorXd& sv, const Eigen::Vector3d& rv_pas) const;

    State initial_rel;
};