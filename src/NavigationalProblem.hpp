#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>

#include "LinearAlgebra.hpp"
#include "GPSHandler.hpp"

class NavigationalProblem {

public:
    NavigationalProblem(const std::string& gnv_filename, const std::string& gps_filename, const GPSHandler& handler);
    void solve(unsigned ti = 0, unsigned tf = 0, bool is_rel = true);
    void errors_norm();
    void errors_prs();
    void errors_rel();

private:
    void load_gnv_data(const std::string& gnv_filename);
    void load_gps_data(const std::string& gps_filename);
    std::vector<double> correct(unsigned prn_id, unsigned gps_time, const std::vector<double>& measurments, bool is_rel);
    std::vector<double> iterative(const std::vector<double>& pseudoranges, const std::vector<std::vector<double>>& gps_positions) const;

    double c = 2.99792458e8;
    double earth_rotation_rate = 7.2921151467e-5;
    double C1 = 2.54572778016;
    double C2 = -1.54572778016;

    double GDOP0 = 5;

    GPSHandler handler;

    std::map<unsigned, std::vector<double>> pos;
    std::map<unsigned, std::vector<double>> vel;

    std::set<unsigned> ts;
    std::map<std::pair<unsigned, unsigned>, std::vector<double>> raw; // [prn_id, gps_time] - [L1_range, L2_range]

    std::map<std::pair<unsigned, unsigned>, double> err_prs;
    std::map<unsigned, std::vector<double>> err_sol;
};
