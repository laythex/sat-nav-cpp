#pragma once

#include <cstdlib>
#include <format>
#include <fstream>

#include <string>

#include "SatNav.hpp"

class PlotterRel {
public:
    PlotterRel(const SatNavRel& sn);

    void plot_errors_norm(double ymin = 0, double ymax = 0);
    void plot_errors_proj(double ymin = 0, double ymax = 0);
    void plot_true_norm(double ymin = 0, double ymax = 0);
    void plot_true_proj(double ymin = 0, double ymax = 0);

    void plot_errors_norm_by_type(IntendedError error, double ymin = 0, double ymax = 0);

private:
    void run_py_plotter(const std::string& arg) const;

    std::vector<double> ECEF_to_geographical(const std::vector<double>& position);

    SatNavRel problem;
    SatNavRel problem_copy;

    char sep = ',';
};
