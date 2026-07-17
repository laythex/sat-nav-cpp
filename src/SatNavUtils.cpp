#include "SatNavUtils.hpp"

double SatNavUtils::ionofree(double l1, double l2) { return l1 * Constants::C1 + l2 * Constants::C2; }

Eigen::Vector3d SatNavUtils::acceleration(const Eigen::Vector3d& rv) {
    double r = rv.norm();
    double z = rv(2);
    Eigen::Vector3d k = Eigen::Vector3d::UnitZ();

    double ir = 1.0 / r;
    double ir2 = ir * ir;
    double ir3 = ir2 * ir;
    double ir5 = ir3 * ir2;

    Eigen::Vector3d a1 = -Constants::mu * ir3 * rv;

    double a2d = -1.5 * Constants::mu * Constants::J2 * Constants::Re * Constants::Re * ir5;
    Eigen::Vector3d a2v = ((1.0 - 5.0 * z * z * ir2) * rv + 2.0 * z * k);
    Eigen::Vector3d a2 = a2d * a2v;

    return a1 + a2;
}

Eigen::Vector3d SatNavUtils::acceleration_rel(const Eigen::Vector3d& Rv, const Eigen::Vector3d rv_pas) {
    return acceleration(rv_pas + Rv) - acceleration(rv_pas);
}

DiscreteState SatNavUtils::step_inc(const DiscreteState& s, double dt, double k) {
    Eigen::Vector3d Av = acceleration(s.position);
    Eigen::Vector3d drv = Av * dt * dt + (Constants::I3 - 2.0 * Constants::Omega * dt) * s.increment * k -
                          Constants::Omega2 * dt * dt * s.position;
    Eigen::Vector3d rv = s.position + drv;
    return {s.time + static_cast<unsigned>(dt), rv, drv / k};
}

DiscreteState SatNavUtils::step_inc_rel(const DiscreteState& s_rel, const Eigen::Vector3d rv_pas, double dt, double k) {
    Eigen::Vector3d Av = acceleration_rel(s_rel.position, rv_pas);
    Eigen::Vector3d dRv = Av * dt * dt + (Constants::I3 - 2.0 * Constants::Omega * dt) * s_rel.increment * k -
                          Constants::Omega2 * dt * dt * s_rel.position;
    Eigen::Vector3d Rv = s_rel.position + dRv;
    return {s_rel.time + static_cast<unsigned>(dt), Rv, dRv / k};
}

State SatNavUtils::step_rk4(const State& s, double dt) {
    Eigen::VectorXd sv(6);
    sv << s.position, s.velocity;

    Eigen::VectorXd k1 = rhs_rk4(sv);
    Eigen::VectorXd k2 = rhs_rk4(sv + k1 * dt / 2.0);
    Eigen::VectorXd k3 = rhs_rk4(sv + k2 * dt / 2.0);
    Eigen::VectorXd k4 = rhs_rk4(sv + k3 * dt);

    sv = sv + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

    return {s.time + static_cast<unsigned>(dt), sv.head(3), sv.tail(3)};
}

State SatNavUtils::step_rk4_rel(const State& s_rel, const State& s_pas, double dt) {
    Eigen::VectorXd sv(6);
    sv << s_rel.position, s_rel.velocity;

    Eigen::Vector3d rv_pas_1 = s_pas.position;
    Eigen::Vector3d rv_pas_2 = step_rk4(s_pas, dt / 2.0).position;
    Eigen::Vector3d rv_pas_4 = step_rk4(s_pas, dt).position;

    Eigen::VectorXd k1 = rhs_rk4_rel(sv, rv_pas_1);
    Eigen::VectorXd k2 = rhs_rk4_rel(sv + k1 * dt / 2.0, rv_pas_2);
    Eigen::VectorXd k3 = rhs_rk4_rel(sv + k2 * dt / 2.0, rv_pas_2);
    Eigen::VectorXd k4 = rhs_rk4_rel(sv + k3 * dt, rv_pas_4);

    sv = sv + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

    return {s_rel.time + static_cast<unsigned>(dt), sv.head(3), sv.tail(3)};
}

Eigen::VectorXd SatNavUtils::rhs_rk4(const Eigen::VectorXd& s) {
    Eigen::Vector3d rv = s.head<3>();
    Eigen::Vector3d vv = s.tail<3>();
    Eigen::Vector3d av = acceleration(rv);

    double ax_ni = 2.0 * Constants::omega * vv[1] + Constants::omega * Constants::omega * rv[0];
    double ay_ni = -2.0 * Constants::omega * vv[0] + Constants::omega * Constants::omega * rv[1];

    Eigen::VectorXd res_sv(6);
    res_sv << vv[0], vv[1], vv[2], av[0] + ax_ni, av[1] + ay_ni, av[2];
    return res_sv;
}

Eigen::VectorXd SatNavUtils::rhs_rk4_rel(const Eigen::VectorXd& s_rel, const Eigen::Vector3d& rv_pas) {
    Eigen::Vector3d Rv = s_rel.head<3>();
    Eigen::Vector3d Vv = s_rel.tail<3>();
    Eigen::Vector3d Av = acceleration_rel(Rv, rv_pas);

    double ax_ni = 2.0 * Constants::omega * Vv[1] + Constants::omega * Constants::omega * Rv[0];
    double ay_ni = -2.0 * Constants::omega * Vv[0] + Constants::omega * Constants::omega * Rv[1];

    Eigen::VectorXd res_sv(6);
    res_sv << Vv[0], Vv[1], Vv[2], Av[0] + ax_ni, Av[1] + ay_ni, Av[2];
    return res_sv;
}

Eigen::Matrix3d SatNavUtils::G(const Eigen::Vector3d& rv) {
    double r = rv.norm();
    double ir2 = 1.0 / (r * r);

    double x = rv(0);
    double y = rv(1);
    double z = rv(2);

    double x2 = x * x;
    double y2 = y * y;
    double z2 = z * z;

    double ac = Constants::mu * ir2 / r;
    double an = 1.5 * Constants::J2 * Constants::Re * Constants::Re * ir2;

    double ft = 5.0 / 3.0;

    double Gxx = -ac * (1.0 - 3.0 * x2 * ir2 + an * (1.0 - 5.0 * (x2 + z2) * ir2 + 35.0 * x2 * z2 * ir2 * ir2));
    double Gyy = -ac * (1.0 - 3.0 * y2 * ir2 + an * (1.0 - 5.0 * (y2 + z2) * ir2 + 35.0 * y2 * z2 * ir2 * ir2));
    double Gzz = -ac * (1.0 - 3.0 * z2 * ir2 + an * (3.0 - 30.0 * z2 * ir2 + 35.0 * z2 * z2 * ir2 * ir2));
    double Gxy = ac * 3.0 * x * y * ir2 * (1.0 + an * ft * (1.0 - 7.0 * z2 * ir2));
    double Gxz = ac * 3.0 * x * z * ir2 * (1.0 + an * ft * (3.0 - 7.0 * z2 * ir2));
    double Gyz = ac * 3.0 * y * z * ir2 * (1.0 + an * ft * (3.0 - 7.0 * z2 * ir2));

    Eigen::Matrix3d G;
    G << Gxx, Gxy, Gxz, Gxy, Gyy, Gyz, Gxz, Gyz, Gzz;
    return G;
}

Eigen::MatrixXd SatNavUtils::iF(const Eigen::Matrix3d& G, double dt) {
    Eigen::Matrix3d iF1 = (Constants::Omega2 - G) * dt * dt;
    Eigen::Matrix3d iF2 = 2.0 * Constants::Omega * dt;
    Eigen::MatrixXd iF(6, 6);
    iF << Constants::I3, -Constants::I3, iF1, Constants::I3 + iF2 - iF1;
    return iF;
}

double SatNavUtils::delta(const Eigen::Vector3d& av, const Eigen::Vector3d& xv, bool use3) {
    double a = av.norm();
    double x = xv.norm();
    double ax = av.dot(xv);

    double ia = 1 / a;
    double ia3 = ia * ia * ia;
    double ia5 = ia3 * ia * ia;

    double d2 = 0.5 * (x * x * ia - ax * ax * ia3);
    double d3 = 0.5 * (ax * ax * ax * ia5 - ax * x * x * ia3);

    return d2 + d3 * use3;
}

Eigen::VectorXd SatNavUtils::shifted_difference(const Eigen::VectorXd& v) {
    if (v.size() <= 1) {
        return Eigen::VectorXd(v.size());
    }

    Eigen::VectorXd shifted(v.size());
    shifted << v.tail(v.size() - 1), v.head(1);
    return v - shifted;
}

Eigen::MatrixX3d SatNavUtils::shifted_difference(const Eigen::MatrixX3d& m) {
    if (m.rows() <= 1) {
        return Eigen::MatrixX3d::Zero(m.rows(), 3);
    }

    Eigen::MatrixX3d shifted(m.rows(), 3);
    shifted << m.bottomRows(m.rows() - 1), m.topRows(1);
    return m - shifted;
}

Eigen::VectorXd SatNavUtils::form_groupped(const std::vector<double>& groupped,
                                           const std::vector<size_t>& present_prns) {
    size_t n = present_prns.size();
    Eigen::VectorXd formed(n);

    for (size_t i = 0; i < n; ++i) {
        size_t prn_index = present_prns[i] - 1;
        Eigen::Index ie = static_cast<Eigen::Index>(i);
        formed(ie) = groupped[prn_index];
    }

    return formed;
}

Eigen::MatrixX3d SatNavUtils::form_groupped(const std::vector<Eigen::Vector3d>& groupped,
                                            const std::vector<size_t>& present_prns) {
    size_t n = present_prns.size();
    Eigen::MatrixX3d formed(n, 3);

    for (size_t i = 0; i < n; ++i) {
        size_t prn_index = present_prns[i] - 1;
        Eigen::Index ie = static_cast<Eigen::Index>(i);
        formed.row(ie) = groupped[prn_index];
    }

    return formed;
}
