#pragma once

#include <iostream>
#include <iomanip>

#include <fstream>
#include <vector>
#include <map>
#include <set>

class NavigationalProblem {

public:
    NavigationalProblem(std::string gnv_filename, std::string gps_filename);
    void solve(unsigned ti, unsigned tf);

private:
    void load_gnv_data(std::string gnv_filename);
    void load_gps_data(std::string gps_filename);

    std::map<unsigned, std::vector<double>> pos;
    std::map<unsigned, std::vector<double>> vel;
    std::set<unsigned> ts;
    std::map<std::pair<unsigned, unsigned>, std::vector<double>> pr;
};