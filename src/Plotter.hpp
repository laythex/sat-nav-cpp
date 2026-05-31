#pragma once

#include <fstream>
#include <cstdlib>

#include <string>
#include <map>

#include "LinAlg.hpp"
#include "SatNav.hpp"

class Plotter {
public:
    Plotter(const SatNav& sn, unsigned ti = 0, unsigned tf = 0);
    void plot_errors_norm(double ymin = 0, double ymax = 0);
    void plot_errors_pr(double ymin = 0, double ymax = 0);
    void plot_gdop(double ymin = 0, double ymax = 0);
    void plot_errors_by_type(IntendedError error = IntendedError::NONE, double ymin = 0, double ymax = 0);

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
    
    unsigned ti, tf;

    char sep = ',';
};
