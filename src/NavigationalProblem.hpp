#pragma once

#include <iostream>
#include <iomanip>

#include <fstream>
#include <vector>
#include <map>
#include <set>

#include "Matrix.hpp"
#include "GPSHandler.hpp"

class NavigationalProblem {

public:
    NavigationalProblem(std::string gnv_filename, std::string gps_filename, const GPSHandler& handler);
    void solve(unsigned ti, unsigned tf);

private:
    void load_gnv_data(std::string gnv_filename);
    void load_gps_data(std::string gps_filename);
    void iterative(const std::vector<double>& PR, const std::map<unsigned, std::vector<double>>& X);

    double c = 2.99792458e8;
    double C1 = 2.54572778016;
    double C2 = -1.54572778016;

    GPSHandler handler;

    std::map<unsigned, std::vector<double>> pos;
    std::map<unsigned, std::vector<double>> vel;
    std::set<unsigned> ts;
    std::map<std::pair<unsigned, unsigned>, std::vector<double>> pr; // [prn_id, gps_time] - [L1_range, L2_range]
};
