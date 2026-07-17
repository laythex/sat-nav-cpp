#include "Config.hpp"
#include "GPSHandler.hpp"
#include "Plotter.hpp"
#include "PlotterRel.hpp"
#include "SatNav.hpp"

int main() {
    Config::ti = 500;
    Config::tf = 86400;

    Config::GDOP0 = 10.0;
    Config::CN0_min_threshold = 25.0;
    Config::CN0_max_threshold = 60.0;
    Config::fadein_threshold = 60;
    Config::hatch_constant = 1 / 10.0;
    Config::solution_tolerance = 0.1;
    Config::residual_threshold = 2.5;
    Config::residual_rel_threshold = 10.0;

    Config::use_kf = true;
    Config::lambda_r = 0.97;
    Config::lambda_v = 0.97;
    Config::lambda_model = 0.98;
    Config::trace_threshold = 30;
    Config::use_rk4 = true;

    Date date = {2005, 12, 10};
    GPSHandler handler = GPSHandler(date);
    SatNav problem1 = SatNav(date, GRACE::A, handler);
    SatNav problem2 = SatNav(date, GRACE::B, handler);

    // Date date = {2023, 3, 23};
    // GPSHandler handler = GPSHandler(date);
    // SatNav problem1 = SatNav(date, GRACE::C, handler);
    // SatNav problem2 = SatNav(date, GRACE::D, handler);

    SatNavRel problem_rel(problem1, problem2);
    PlotterRel plotter(problem_rel);

    double m = 3.0;
    plotter.plot_errors_norm(0, m);
    // plotter.plot_errors_norm_by_type(IntendedError::DELTA3, 0, 10.0);
    // plotter.plot_errors_proj(-m, m);
    // plotter.plot_true_norm(0.0, 200.0);

    // Date date = {2026, 3, 1};
    // GPSHandler handler = GPSHandler(date);
    // SatNav problem1 = SatNav(date, SWARM::C, handler);
    // SatNav problem2 = SatNav(date, SWARM::A, handler);

    // SatNavRel problem_rel(problem1, problem2);
    // PlotterRel plotter(problem_rel, 1000, 4000);

    // double m = 0;
    // plotter.plot_errors_norm(0, m);
    // plotter.plot_errors_proj(-m, m);
    // plotter.plot_true_norm();
    // plotter.plot_true_proj();

    return 0;
}