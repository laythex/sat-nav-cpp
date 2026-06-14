#include "Propagator.hpp"

Propagator::Propagator(const State& initial, bool useJ2) : initial(initial), useJ2(useJ2) {}

State Propagator::step_rk4(double step) {
    Eigen::VectorXd sv(6);
    sv << initial.position[0], initial.position[1], initial.position[2], initial.velocity[0], initial.velocity[1],
        initial.velocity[2];

    Eigen::VectorXd k1 = rhs_rk4(sv);
    Eigen::VectorXd k2 = rhs_rk4(sv + k1 * step / 2.0);
    Eigen::VectorXd k3 = rhs_rk4(sv + k2 * step / 2.0);
    Eigen::VectorXd k4 = rhs_rk4(sv + k3 * step);

    sv = sv + step / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

    State s = {initial.time + static_cast<unsigned>(step), {sv[0], sv[1], sv[2]}, {sv[3], sv[4], sv[5]}};
    initial = s;
    return s;
}

State Propagator::step_inc(double step) {
    Eigen::Vector3d rv = {initial.position[0], initial.position[1], initial.position[2]};
    Eigen::Vector3d dv = {initial.velocity[0], initial.velocity[1], initial.velocity[2]};

    Eigen::Vector3d acc = acceleration(rv);

    dv = acc * step * step + (Constants::I3 - 2.0 * Constants::Omega * step) * dv -
         Constants::Omega * Constants::Omega * step * step * rv;
    rv = rv + dv;

    State s = {initial.time + static_cast<unsigned>(step), {rv[0], rv[1], rv[2]}, {dv[0], dv[1], dv[2]}};
    initial = s;
    return s;
}

Eigen::Vector3d Propagator::acceleration(const Eigen::Vector3d& rv) const {
    double r = rv.norm();
    double z = rv[2];
    Eigen::Vector3d k = {0.0, 0.0, 1.0};

    Eigen::Vector3d a1 = -Constants::mu / (r * r * r) * rv;

    double a2d = -3.0 * Constants::mu * Constants::J2 * Constants::Re * Constants::Re / (2.0 * r * r * r * r * r);
    Eigen::Vector3d a2v = ((1.0 - 5.0 * z * z / (r * r)) * rv + 2.0 * z * k);
    Eigen::Vector3d a2 = a2d * a2v;

    if (!useJ2) {
        return a1;
    }

    return a1 + a2;
}

Eigen::VectorXd Propagator::rhs_rk4(const Eigen::VectorXd& sv) const {
    Eigen::Vector3d rv = {sv[0], sv[1], sv[2]};
    Eigen::Vector3d acc = acceleration(rv);

    double ax_ni = 2.0 * Constants::omega * sv[4] + Constants::omega * Constants::omega * sv[0];
    double ay_ni = -2.0 * Constants::omega * sv[3] + Constants::omega * Constants::omega * sv[1];

    Eigen::VectorXd res_sv(6);
    res_sv << sv[3], sv[4], sv[5], acc[0] + ax_ni, acc[1] + ay_ni, acc[2];
    return res_sv;
}
