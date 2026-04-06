#pragma once

#include <fstream>
#include <cstdlib>

#include <string>
#include <map>

#include "LinAlg.hpp"
#include "SatNav.hpp"

class Plotter {
public:
    Plotter(const SatNav& sn, unsigned time_initial = 0, unsigned time_final = 0);
    void plot_errors_norm(double ymin = 0, double ymax = 0);
    void plot_errors_pr(double ymin = 0, double ymax = 0);
    void plot_gdop(double ymin = 0, double ymax = 0);
    void plot_errors_by_type(char error_type, double ymin = 0, double ymax = 0);

    void plot_map_norm(const std::string& mode, double threshold = 0);
    void plot_map_iono(const std::string& mode, double threshold = 0);

    void plot_acc_lin_proj(double ymin = 0, double ymax = 0);
    void plot_acc_ang_proj(double ymin = 0, double ymax = 0);

private:
    void run_py_plotter(const std::string& arg) const;
    void run_py_mapper(const std::string& arg) const;

    std::vector<double> ECEF_to_geographycal(const std::vector<double>& position); // сделать отдельную структуру гео координат?

    SatNav problem;
    SatNav problem_copy;

    char sep = ',';
    
    unsigned time_initial_, time_final_;

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
