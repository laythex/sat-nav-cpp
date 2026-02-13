#pragma once

#include <fstream>
#include <cstdlib>

#include <string>

#include "LinAlg.hpp"
#include "SatNav.hpp"

// Объединить в один класс через темплейт?

class PlotterRel {
public:
    PlotterRel(const SatNavRel& sn, double ti = 0, double tf = 0);
    void plot_errors_norm(double ymin = 0, double ymax = 0);

private:
    void run_py_plotter(const std::string& arg) const;

    std::vector<double> ECEF_to_geographical(const std::vector<double>& position);

    SatNavRel problem;

    char sep = ',';
    
    unsigned ti, tf;
};
