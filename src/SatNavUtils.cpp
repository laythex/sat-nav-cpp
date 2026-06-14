#include "SatNavUtils.hpp"

bool SatNavUtils::check_time(unsigned t, unsigned ti, unsigned tf, unsigned t0) {
    if (ti == 0 && tf == 0) {
        return true;
    }

    unsigned td = t - t0;
    if (td >= ti && td < tf) {
        return true;
    }

    return false;
}

double SatNavUtils::ionofree(double l1, double l2) { return l1 * Constants::C1 + l2 * Constants::C2; }

double SatNavUtils::calculate_delta(const Eigen::Vector3d& av, const Eigen::Vector3d& xv) {
    double a = av.norm();
    double x = xv.norm();
    double ax = av.dot(xv);

    double ia = 1 / a;
    double ia3 = ia * ia * ia;
    double ia5 = ia3 * ia * ia;

    double d2 = 0.5 * (x * x * ia - ax * ax * ia3);
    double d3 = 0.5 * (ax * ax * ax * ia5 - ax * x * x * ia3);

    return d2 + d3;
}

Eigen::VectorXd SatNavUtils::shifted_difference(const Eigen::VectorXd& v) {
    Eigen::VectorXd shifted(v.size());
    shifted << v.tail(v.size() - 1), v.head(1);
    return v - shifted;
}

Eigen::MatrixX3d SatNavUtils::shifted_difference(const Eigen::MatrixX3d& m) {
    Eigen::MatrixX3d shifted(m.rows(), 3);
    shifted << m.bottomRows(m.rows() - 1), m.topRows(1);
    return m - shifted;
}