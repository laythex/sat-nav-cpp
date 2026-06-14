#pragma once

#include <cmath>

#include <map>
#include <string>
#include <vector>

#include <iostream>

#include "DataParser.hpp"

class GPSHandler {

public:
    GPSHandler(const Date& date);

    double get_clock_error(unsigned prn_id, double gps_time) const;
    double get_relativistic_error(unsigned prn_id, double gps_time) const;
    State get_state(unsigned prn_id, double gps_time) const;

private:
    const Ephemeris& select_ephemeris(unsigned prn_id, double t_sv) const;

    static double gps_to_sv(double gps_time);

    double mu_sqrt = 1.99649818432e7;
    double Omega_e_dot = 7.2921151467e-5;
    double F = -4.442807633e-10;

    std::vector<std::vector<Ephemeris>> ephemeris;
};