#pragma once

#include <iostream>

#include <fstream>
#include <iomanip>

#include <string>
#include <vector>
#include <map>
#include <set>

struct Ephemeris {
    double a_f0, a_f1, a_f2;
    double M_0, delta_n, e, A_sqrt, Omega_0, i_0, omega, Omega_dot, IDOT;
    double C_uc, C_us, C_rc, C_rs, C_ic, C_is;
};

namespace DataParser {
    void load_fo_gnv_data(const std::string& gnv_filename, 
                          std::map<unsigned, std::vector<double>>& pos, 
                          std::map<unsigned, std::vector<double>>& vel);

    void load_fo_gps_data(const std::string& gps_filename, 
                          std::set<unsigned>& ts_raw,
                          std::map<std::pair<unsigned, unsigned>, std::vector<double>>& raw);
    
    void load_gnv_data(const std::string& gnv_filename, 
                       std::map<unsigned, std::vector<double>>& pos, 
                       std::map<unsigned, std::vector<double>>& vel);
    
    void load_gps_data(const std::string& gps_filename, 
                       std::set<unsigned>& ts_raw,
                       std::map<std::pair<unsigned, unsigned>, std::vector<double>>& raw);

    void load_brdc_data(const std::string& brdc_filename, 
                        std::map<unsigned, std::vector<unsigned>>& ts_ephs, 
                        std::map<std::pair<unsigned, unsigned>, Ephemeris>& ephs);
};