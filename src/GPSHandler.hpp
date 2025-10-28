#pragma once

#include <iostream>
#include <iomanip>

#include <cmath>
#include <fstream>
#include <vector>
#include <map>

class GPSHandler {

public:
    GPSHandler(std::string rnx_filename);

    std::vector<double> interp_state(unsigned prn_id, double gps_time);

    static unsigned convert_id(unsigned prn_id);

private:
    void load_rnx_data(std::string rnx_filename);

    unsigned get_ephemeris(unsigned prn_id, double t_sv);

    double c = 2.99792458e8;
    double mu = 3.986005e14;
    double Omega_e_dot = 7.2921151467e-5;
    double F = -4.442807633e-10;

    std::map<unsigned, std::vector<unsigned>> ephemeris;

    std::map<std::pair<unsigned, unsigned>, double> map_a_f0;
    std::map<std::pair<unsigned, unsigned>, double> map_a_f1;
    std::map<std::pair<unsigned, unsigned>, double> map_a_f2;
    std::map<std::pair<unsigned, unsigned>, double> map_M_0;
    std::map<std::pair<unsigned, unsigned>, double> map_delta_n;
    std::map<std::pair<unsigned, unsigned>, double> map_e;
    std::map<std::pair<unsigned, unsigned>, double> map_A_sqrt;
    std::map<std::pair<unsigned, unsigned>, double> map_Omega_0;
    std::map<std::pair<unsigned, unsigned>, double> map_i_0;
    std::map<std::pair<unsigned, unsigned>, double> map_omega;
    std::map<std::pair<unsigned, unsigned>, double> map_Omega_dot;
    std::map<std::pair<unsigned, unsigned>, double> map_IDOT;
    std::map<std::pair<unsigned, unsigned>, double> map_C_uc;
    std::map<std::pair<unsigned, unsigned>, double> map_C_us;
    std::map<std::pair<unsigned, unsigned>, double> map_C_rc;
    std::map<std::pair<unsigned, unsigned>, double> map_C_rs;
    std::map<std::pair<unsigned, unsigned>, double> map_C_ic;
    std::map<std::pair<unsigned, unsigned>, double> map_C_is;
};