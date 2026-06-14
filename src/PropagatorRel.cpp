#include "PropagatorRel.hpp"

PropagatorRel::PropagatorRel(const State& initial_rel) : initial_rel(initial_rel) {}

State PropagatorRel::step_rk4_rel(double step, const State& s_pas) {
    Eigen::VectorXd sv(6);
    sv << initial_rel.position[0], initial_rel.position[1], initial_rel.position[2], initial_rel.velocity[0],
        initial_rel.velocity[1], initial_rel.velocity[2];

    Propagator propagator_pas(s_pas);
    Eigen::Vector3d rv_pas_1 = Eigen::Vector3d::Map(s_pas.position.data());
    Eigen::Vector3d rv_pas_2 = Eigen::Vector3d::Map(propagator_pas.step_rk4(step / 2.0).position.data());
    Eigen::Vector3d rv_pas_4 = Eigen::Vector3d::Map(propagator_pas.step_rk4(step / 2.0).position.data());

    Eigen::VectorXd k1 = rhs_rk4_rel(sv, rv_pas_1);
    Eigen::VectorXd k2 = rhs_rk4_rel(sv + k1 * step / 2.0, rv_pas_2);
    Eigen::VectorXd k3 = rhs_rk4_rel(sv + k2 * step / 2.0, rv_pas_2);
    Eigen::VectorXd k4 = rhs_rk4_rel(sv + k3 * step, rv_pas_4);

    sv = sv + step / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

    State s_rel = {initial_rel.time + static_cast<unsigned>(step), {sv[0], sv[1], sv[2]}, {sv[3], sv[4], sv[5]}};
    initial_rel = s_rel;
    return s_rel;
}

State PropagatorRel::step_inc_rel(double step, const State& s_pas) {
    Eigen::Vector3d rv = {initial_rel.position[0], initial_rel.position[1], initial_rel.position[2]};
    Eigen::Vector3d dv = {initial_rel.velocity[0], initial_rel.velocity[1], initial_rel.velocity[2]};

    Eigen::Vector3d rv_pas = Eigen::Vector3d::Map(s_pas.position.data());
    Eigen::Vector3d acc = acceleration_rel(rv, rv_pas);

    dv = acc * step * step + (Constants::I3 - 2.0 * Constants::Omega * step) * dv -
         Constants::Omega * Constants::Omega * step * step * rv;
    rv = rv + dv;

    State s_rel = {initial_rel.time + static_cast<unsigned>(step), {rv[0], rv[1], rv[2]}, {dv[0], dv[1], dv[2]}};
    initial_rel = s_rel;
    return s_rel;
}

Eigen::Vector3d PropagatorRel::acceleration(const Eigen::Vector3d& rv) const {
    double r = rv.norm();
    double z = rv[2];
    Eigen::Vector3d k = {0.0, 0.0, 1.0};

    Eigen::Vector3d a1 = -Constants::mu / (r * r * r) * rv;

    double a2d = -3.0 * Constants::mu * Constants::J2 * Constants::Re * Constants::Re / (2.0 * r * r * r * r * r);
    Eigen::Vector3d a2v = ((1.0 - 5.0 * z * z / (r * r)) * rv + 2.0 * z * k);
    Eigen::Vector3d a2 = a2d * a2v;

    return a1 + a2;
}

Eigen::Vector3d PropagatorRel::acceleration_rel(const Eigen::Vector3d& Rv, const Eigen::Vector3d& rv_pas) const {
    return acceleration(rv_pas + Rv) - acceleration(rv_pas);
}

Eigen::VectorXd PropagatorRel::rhs_rk4_rel(const Eigen::VectorXd& sv, const Eigen::Vector3d& rv_pas) const {
    Eigen::Vector3d Rv = sv.head<3>();
    Eigen::Vector3d Vv = sv.tail<3>();
    Eigen::Vector3d Av = acceleration_rel(Rv, rv_pas);

    double ax_ni = 2.0 * Constants::omega * Vv[1] + Constants::omega * Constants::omega * Rv[0];
    double ay_ni = -2.0 * Constants::omega * Vv[0] + Constants::omega * Constants::omega * Rv[1];

    Eigen::VectorXd res_sv(6);
    res_sv << Vv[0], Vv[1], Vv[2], Av[0] + ax_ni, Av[1] + ay_ni, Av[2];
    return res_sv;
}
