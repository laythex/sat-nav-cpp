#pragma once

#include <cmath>

#include <string>
#include <vector>
#include <map>

#include <iostream>

#include "DataParser.hpp"

class GPSHandler {

public:
    GPSHandler(const std::string& brdc_filename);
    double get_clock_error(unsigned prn_id, double grace_time);
    GPSState get_state(unsigned prn_id, double grace_time);

private:
    double grace_to_sv(double grace_time) const;
    unsigned select_ephemeris(unsigned prn_id, double t_sv);

    double mu_sqrt = 1.99649818432e7;
    double Omega_e_dot = 7.2921151467e-5;
    double F = -4.442807633e-10;

    std::map<unsigned, std::vector<unsigned>> ts_ephs; // prn_id - [eph_t1, eph_t2, ...]
    std::map<std::pair<unsigned, unsigned>, Ephemeris> ephs; // [prn_id, eph_t] - Ephemeris
};