#pragma once

#include <cstdlib>
#include <format>
#include <fstream>

#include <map>
#include <string>

#include "SatNav.hpp"

class Plotter {
public:
    Plotter(const SatNav& sn, unsigned ti = 0, unsigned tf = 0);
    void plot_errors_norm(double ymin = 0, double ymax = 0);
    void plot_errors_pr(double ymin = 0, double ymax = 0);
    void plot_gdop(double ymin = 0, double ymax = 0);
    void plot_errors_norm_by_type(IntendedError error = IntendedError::NONE, double ymin = 0, double ymax = 0);
    void plot_errors_proj_by_type(IntendedError error = IntendedError::NONE, double ymin = 0, double ymax = 0);
    void plot_errors_pr_by_type(IntendedError error = IntendedError::NONE, double ymin = 0, double ymax = 0);

    void plot_map_iono(const std::string& mode, double threshold = 0);

    void plot_acc_lin_proj(double ymin = 0, double ymax = 0);
    void plot_acc_ang_proj(double ymin = 0, double ymax = 0);

private:
    void run_py_plotter(const std::string& arg) const;
    void run_py_mapper(const std::string& arg) const;

    void reset_problem_copy();

    Eigen::Vector3d ECEF_to_geographycal(const Eigen::Vector3d& position);

    SatNav problem;
    SatNav problem_copy;

    unsigned ti, tf;

    char sep = ',';
};
