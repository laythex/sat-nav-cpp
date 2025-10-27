#include "GPSHandler.hpp"

GPSHandler::GPSHandler(string rnx_data_filename) {

}

GPSHandler::load_rnx_data(string rnx_data_filename) {

}

std::vector<double> GPSHandler::interp_state(unsigned prn_id, unsigned gps_time) {

    unsigned timestamp = 0;
    std::pair key = std::make_pair(prn_id, timestamp);

    double M_0 = map_M_0[key];
    double delta_n = map_delta_n[key];
    double e = map_e[key];
    double A_sqrt = map_A_sqrt[key];
    double Omega_0 = map_Omega_0[key];
    double i_0 = map_i_0[key];
    double omega = map_omega[key];
    double Omega_dot = map_Omega_dot[key];
    double IDOT = map_IDOT[key];
    double C_uc = map_C_uc[key];
    double C_us = map_C_us[key];
    double C_rc = map_C_rc[key];
    double C_rs = map_C_rs[key];
    double C_ic = map_C_ic[key];
    double C_is = map_C_is[key];
    double t_oe = map_t_oe[key];

    double a_f0 = map_a_f0[key];
    double a_f1 = map_a_f1[key];
    double a_f2 = map_a_f2[key];

    double t = t_sv - t_oe;
    double delta_t_sv = a_f0 + a_f1 * t + a_f2 * t * t;

    double A = A_sqrt * A_sqrt;
    double n_0 = sqrt(mu) / (A_sqrt * A_sqrt * A_sqrt);
    double n = n_0 + delta_n;
    double M = M_0 + n * t;

    double E = M;
    for (unsigned i = 0; i < 3; i++) {
        E += (M - E + e * sin(E)) / (1 - e * cos(E));
    }

    double nu = 2 * arctan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
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

    double Omega_dot = Omega_dot - Omega_e_dot
    
    double x_prime_dot = r_dot * cos(u) - r * u_dot * sin(u);
    double y_prime_dot = r_dot * sin(u) + r * u_dot * cos(u);

    double vx = -x_prime * Omega_dot * sin(Omega) + x_prime_dot * cos(Omega) - y_prime_dot * sin(Omega) * cos(i) - 
                y_prime * (Omega_dot * cos(Omega) * cos(i) - i_dot * sin(Omega) * sin(i));
    double vy = x_prime * Omega_dot * cos(Omega) + x_prime_dot * sin(Omega) + y_prime_dot * cos(Omega) * cos(i) - 
                y_prime * (Omega_dot * sin(Omega) * cos(i) + i_dot * cos(Omega) * sin(i));
    double vz = y_prime_dot * sin(i) + y_prime * i_dot * cos(i);

    double delta_t_r = F * e * A_sqrt * sin(E);
    t += delta_t_r;

    return {t, x, y, z, vx, vy, vz};
}

