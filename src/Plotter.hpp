#pragma once

#include <fstream>
#include <cstdlib>

#include <string>
#include <map>

#include "LinAlg.hpp"
#include "SatNav.hpp"

class Plotter {
public:
    Plotter(const SatNav& sn);
    void plot_errors_norm(double ymin = 0, double ymax = 0);
    void plot_errors_pr(double ymin = 0, double ymax = 0);
    void plot_errors_by_type(char error_type, double ymin = 0, double ymax = 0);

    void plot_map_iono();

private:
    void run_py_plotter(const std::string& arg) const;
    void run_py_map_plotter(const std::string& arg) const;

    std::vector<double> ECEF_to_geographycal(const std::vector<double>& position); // сделать отдельную структура гео координат

    SatNav problem;
    SatNav problem_copy;

    char sep = ',';

    std::map<char, std::string> error_names = {
        {'F', "fading"},
        {'G', "gdop"},
        {'H', "hatch"},
        {'I', "iono"},
        {'L', "low"},
        {'R', "relativistic"},
        {'s', "sagnac"},
        {'S', "snr"}
    };
    std::map<char, std::string> error_titles = {
        {'G', "gdop"},
        {'H', "hatch"},
        {'I', "iono"},
        {'L', "low"},
        {'R', "relativistic"},
        {'S', "snr"}
    };
};
