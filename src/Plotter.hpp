#pragma once

#include <fstream>
#include <cstdlib>

#include <string>

#include "LinAlg.hpp"
#include "SatNav.hpp"

class Plotter {
public:
    Plotter(SatNav& problem);
    void plot_errors_norm(double ymin = 0, double ymax = 0);
    void plot_errors_pr(double ymin = 0, double ymax = 0);
    void plot_errors_by_type(char error_type, double ymin = 0, double ymax = 0);

    void plot_map_iono();

private:
    void run_py_plotter(const std::string& arg) const;


    SatNav& problem;

    char sep = ',';
};
