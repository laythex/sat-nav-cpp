#include "GPSHandler.hpp"

GPSHandler::GPSHandler(const std::string& brdc_filename) {
    DataParser::load_brdc_data(brdc_filename, ts_ephs, ephs);
}

double GPSHandler::grace_to_sv(double grace_time) const {
    unsigned week_number = (grace_time + 630763200) / 604800;
    unsigned week_start = week_number * 604800 - 630763200;
    return grace_time - week_start;
}

unsigned GPSHandler::select_ephemeris(unsigned prn_id, double t_sv) {
    unsigned n = ts_ephs[prn_id].size();

    unsigned lo = 0, hi = n, mid;
    while (lo < hi) {
        mid = lo + (hi - lo) / 2;
        if (ts_ephs[prn_id][mid] < t_sv) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    if (lo == 0) return ts_ephs[prn_id][0];
    if (lo == n) return ts_ephs[prn_id].back();

    unsigned left = ts_ephs[prn_id][lo - 1];
    unsigned right = ts_ephs[prn_id][lo];
    return (t_sv - left < right - t_sv) ? left : right;
}

double GPSHandler::get_clock_error(unsigned prn_id, double grace_time) {
    double t_sv = grace_to_sv(grace_time);
    unsigned t_oe = select_ephemeris(prn_id, t_sv);
    std::pair key = {prn_id, t_oe};

    double a_f0 = ephs[key].a_f0;
    double a_f1 = ephs[key].a_f1;
    double a_f2 = ephs[key].a_f2;

    double t = t_sv - t_oe;
    double delta_t_sv = a_f0 + a_f1 * t + a_f2 * t * t;

    return delta_t_sv;
}

GPSState GPSHandler::get_state(unsigned prn_id, double grace_time) {
    double t_sv = grace_to_sv(grace_time); // а если не грейс?
    unsigned t_oe = select_ephemeris(prn_id, t_sv);
    std::pair key = {prn_id, t_oe};

    double M_0 = ephs[key].M_0;
    double delta_n = ephs[key].delta_n;
    double e = ephs[key].e;
    double A_sqrt = ephs[key].A_sqrt;
    double Omega_0 = ephs[key].Omega_0;
    double i_0 = ephs[key].i_0;
    double omega = ephs[key].omega;
    double Omega_dot = ephs[key].Omega_dot;
    double IDOT = ephs[key].IDOT;
    double C_uc = ephs[key].C_uc;
    double C_us = ephs[key].C_us;
    double C_rc = ephs[key].C_rc;
    double C_rs = ephs[key].C_rs;
    double C_ic = ephs[key].C_ic;
    double C_is = ephs[key].C_is;

    double t = t_sv - t_oe;

    double A = A_sqrt * A_sqrt;
    double n_0 = mu_sqrt / (A_sqrt * A_sqrt * A_sqrt);
    double n = n_0 + delta_n;
    double M = M_0 + n * t;

    double E = M;
    for (unsigned i = 0; i < 3; i++) {
        E += (M - E + e * sin(E)) / (1 - e * cos(E));
    }

    double nu = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
    double Phi = nu + omega;

    double delta_u = C_us * sin(2 * Phi) + C_uc * cos(2 * Phi);
    double delta_r = C_rs * sin(2 * Phi) + C_rc * cos(2 * Phi);
    double delta_i = C_is * sin(2 * Phi) + C_ic * cos(2 * Phi);

    double u = Phi + delta_u;
    double r = A * (1 - e * cos(E)) + delta_r;
    double i = i_0 + delta_i + IDOT * t;

    double x_prime = r * cos(u);
    double y_prime = r * sin(u);
    double Omega = Omega_0 + (Omega_dot - Omega_e_dot) * t - Omega_e_dot * t_oe;

    double x = x_prime * cos(Omega) - y_prime * cos(i) * sin(Omega);
    double y = x_prime * sin(Omega) + y_prime * cos(i) * cos(Omega);
    double z = y_prime * sin(i);

    double E_dot = n / (1 - e * cos(E));
    double nu_dot = E_dot * sqrt(1 - e * e) / (1 - e * cos(E));
    double i_dot = IDOT + 2 * nu_dot * (C_is * cos(2 * Phi) - C_ic * sin(2 * Phi));
    double u_dot = nu_dot + 2 * nu_dot * (C_us * cos(2 * Phi) - C_uc * sin(2 * Phi));
    double r_dot = e * A * E_dot * sin(E) + 2 * nu_dot * (C_rs * cos(2 * Phi) - C_rc * sin(2 * Phi));

    double Omega_k_dot = Omega_dot - Omega_e_dot;
    
    double x_prime_dot = r_dot * cos(u) - r * u_dot * sin(u);
    double y_prime_dot = r_dot * sin(u) + r * u_dot * cos(u);

    double vx = -x_prime * Omega_k_dot * sin(Omega) + x_prime_dot * cos(Omega) - y_prime_dot * sin(Omega) * cos(i) - 
                y_prime * (Omega_k_dot * cos(Omega) * cos(i) - i_dot * sin(Omega) * sin(i));
    double vy = x_prime * Omega_k_dot * cos(Omega) + x_prime_dot * sin(Omega) + y_prime_dot * cos(Omega) * cos(i) - 
                y_prime * (Omega_k_dot * sin(Omega) * cos(i) + i_dot * cos(Omega) * sin(i));
    double vz = y_prime_dot * sin(i) + y_prime * i_dot * cos(i);

    double delta_t_r = F * e * A_sqrt * sin(E);

    return {{0, {x, y, z}, {vx, vy, vz}}, delta_t_r};
}

