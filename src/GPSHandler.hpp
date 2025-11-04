#pragma once

#include <cmath>
#include <fstream>
#include <vector>
#include <map>

#include <iostream>

struct Ephemeris {
    double a_f0, a_f1, a_f2;
    double M_0, delta_n, e, A_sqrt, Omega_0, i_0, omega, Omega_dot, IDOT;
    double C_uc, C_us, C_rc, C_rs, C_ic, C_is;
};


class GPSHandler {

public:
    GPSHandler(std::string rnx_filename);
    double get_clock_error(unsigned prn_id, double gps_time);
    std::vector<double> get_state(unsigned prn_id, double gps_time);

private:
    void load_rnx_data(std::string rnx_filename);
    unsigned select_ephemeris(unsigned prn_id, double t_sv);

    unsigned conversion = 630763200;
    double mu_sqrt = 1.99649818432e7;
    double Omega_e_dot = 7.2921151467e-5;
    double F = -4.442807633e-10;

    std::map<unsigned, std::vector<unsigned>> ts_ephs; // prn_id - [eph_t1, eph_t2, ...]
    std::map<std::pair<unsigned, unsigned>, Ephemeris> ephs; // [prn_id, eph_t] - Ephemeris
};