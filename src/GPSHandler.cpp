#include "GPSHandler.hpp"

GPSHandler::GPSHandler(const Date& date) { ephemeris = DataParser::load_brdc_data(date); }

double GPSHandler::get_clock_error(unsigned prn_id, double time) const {
    double t_sv = gps_to_sv(time);
    const Ephemeris& eph = select_ephemeris(prn_id, t_sv);

    double t = t_sv - eph.t_oe;
    return eph.a_f0 + eph.a_f1 * t + eph.a_f2 * t * t;
}

double GPSHandler::get_relativistic_error(unsigned prn_id, double time) const {
    double t_sv = gps_to_sv(time);
    const Ephemeris& eph = select_ephemeris(prn_id, t_sv);

    double t = t_sv - eph.t_oe;

    double n_0 = mu_sqrt / (eph.A_sqrt * eph.A_sqrt * eph.A_sqrt);
    double n = n_0 + eph.delta_n;
    double M = eph.M_0 + n * t;

    double E = M;
    for (unsigned i = 0; i < 3; ++i) {
        E += (M - E + eph.e * sin(E)) / (1 - eph.e * cos(E));
    }

    return -F * eph.e * eph.A_sqrt * sin(E);
}

State GPSHandler::get_state(unsigned prn_id, double time) const {
    double t_sv = gps_to_sv(time);
    const Ephemeris& eph = select_ephemeris(prn_id, t_sv);

    double t = t_sv - eph.t_oe;

    double A = eph.A_sqrt * eph.A_sqrt;
    double n_0 = mu_sqrt / (eph.A_sqrt * eph.A_sqrt * eph.A_sqrt);
    double n = n_0 + eph.delta_n;
    double M = eph.M_0 + n * t;

    double E = M;
    for (unsigned i = 0; i < 3; ++i) {
        E += (M - E + eph.e * sin(E)) / (1 - eph.e * cos(E));
    }

    double nu = 2 * atan(sqrt((1 + eph.e) / (1 - eph.e)) * tan(E / 2));
    double Phi = nu + eph.omega;

    double delta_u = eph.C_us * sin(2 * Phi) + eph.C_uc * cos(2 * Phi);
    double delta_r = eph.C_rs * sin(2 * Phi) + eph.C_rc * cos(2 * Phi);
    double delta_i = eph.C_is * sin(2 * Phi) + eph.C_ic * cos(2 * Phi);

    double u = Phi + delta_u;
    double r = A * (1 - eph.e * cos(E)) + delta_r;
    double i = eph.i_0 + delta_i + eph.IDOT * t;

    double x_prime = r * cos(u);
    double y_prime = r * sin(u);
    double Omega = eph.Omega_0 + (eph.Omega_dot - Omega_e_dot) * t - Omega_e_dot * eph.t_oe;

    double x = x_prime * cos(Omega) - y_prime * cos(i) * sin(Omega);
    double y = x_prime * sin(Omega) + y_prime * cos(i) * cos(Omega);
    double z = y_prime * sin(i);

    double E_dot = n / (1 - eph.e * cos(E));
    double nu_dot = E_dot * sqrt(1 - eph.e * eph.e) / (1 - eph.e * cos(E));
    double i_dot = eph.IDOT + 2 * nu_dot * (eph.C_is * cos(2 * Phi) - eph.C_ic * sin(2 * Phi));
    double u_dot = nu_dot + 2 * nu_dot * (eph.C_us * cos(2 * Phi) - eph.C_uc * sin(2 * Phi));
    double r_dot = eph.e * A * E_dot * sin(E) + 2 * nu_dot * (eph.C_rs * cos(2 * Phi) - eph.C_rc * sin(2 * Phi));

    double Omega_k_dot = eph.Omega_dot - Omega_e_dot;

    double x_prime_dot = r_dot * cos(u) - r * u_dot * sin(u);
    double y_prime_dot = r_dot * sin(u) + r * u_dot * cos(u);

    double vx = -x_prime * Omega_k_dot * sin(Omega) + x_prime_dot * cos(Omega) - y_prime_dot * sin(Omega) * cos(i) -
                y_prime * (Omega_k_dot * cos(Omega) * cos(i) - i_dot * sin(Omega) * sin(i));
    double vy = x_prime * Omega_k_dot * cos(Omega) + x_prime_dot * sin(Omega) + y_prime_dot * cos(Omega) * cos(i) -
                y_prime * (Omega_k_dot * sin(Omega) * cos(i) + i_dot * cos(Omega) * sin(i));
    double vz = y_prime_dot * sin(i) + y_prime * i_dot * cos(i);

    return {0, {x, y, z}, {vx, vy, vz}};
}

const Ephemeris& GPSHandler::select_ephemeris(unsigned prn_id, double t_sv) const {
    unsigned prn_index = prn_id - 1;
    size_t n = ephemeris[prn_index].size();

    size_t lo = 0, hi = n, mid;
    while (lo < hi) {
        mid = lo + (hi - lo) / 2;
        if (ephemeris[prn_index][mid].t_oe < t_sv) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    if (lo == 0) return ephemeris[prn_index][0];
    if (lo == n) return ephemeris[prn_index].back();

    const Ephemeris& left = ephemeris[prn_index][lo - 1];
    const Ephemeris& right = ephemeris[prn_index][lo];
    return (t_sv - left.t_oe < right.t_oe - t_sv) ? left : right;
}

double GPSHandler::gps_to_sv(double gps_time) {
    unsigned week_number = static_cast<unsigned>(gps_time / 604800);
    return gps_time - week_number * 604800;
}
