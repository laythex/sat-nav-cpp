#include "Config.hpp"
#include "GPSHandler.hpp"
#include "Plotter.hpp"
#include "PlotterRel.hpp"
#include "SatNav.hpp"

int main() {
    /*
    TODO:
    все string через format
    Raw Measurement хранит C1C, L1C и т.д.
    темплейт бинарного поиска
    */

    Config::ti = 0;
    Config::tf = 86400;

    Config::GDOP0 = 10.0;
    Config::CN0_min_threshold = 32.0;
    Config::CN0_max_threshold = 60.0;
    Config::fadein_threshold = 60;
    Config::hatch_constant = 1 / 10.0;
    Config::solution_tolerance = 0.1;

    Config::use_kf = false;

    // Date date = {2005, 12, 10}; GRACE sat = GRACE::A;
    // Date date = {2005, 12, 10}; GRACE sat = GRACE::B;
    // Date date = {2019, 10, 21};
    // Date date = {2020, 5, 20}; GRACE sat = GRACE::C;
    // Date date = {2021, 3, 30}; GRACE sat = GRACE::C;
    // Date date = {2022, 9, 18}; GRACE sat = GRACE::C;
    // Date date = {2023, 3, 23};
    // Date date = {2024, 4, 3}; GRACE sat = GRACE::C;

    // GRACE sat = GRACE::C;

    // GPSHandler handler = GPSHandler(date);
    // SatNav problem = SatNav(date, sat, handler);

    // double m = 10;
    // Plotter plotter(problem);
    // plotter.plot_errors_pr(-m, m);
    // plotter.plot_errors_norm(0, m);

    // double n = 40;
    // plotter.plot_errors_proj_by_type(IntendedError::RELATIVISTIC, -n, n);
    // plotter.plot_errors_pr_by_type(IntendedError::RELATIVISTIC, -n, n);
    // plotter.plot_errors_proj_by_type(IntendedError::IONOSPHERIC, -n, n);
    // plotter.plot_errors_pr_by_type(IntendedError::IONOSPHERIC, -n, n);
    // plotter.plot_errors_proj_by_type(IntendedError::SAGNAC, -n, n);
    // plotter.plot_errors_pr_by_type(IntendedError::SAGNAC, -n, n);
    // plotter.plot_errors_proj_by_type(IntendedError::NOISY, -n, n);
    // plotter.plot_errors_pr_by_type(IntendedError::NOISY, -n, n);
    // plotter.plot_errors_proj_by_type(IntendedError::HATCH, -n, n);
    // plotter.plot_errors_pr_by_type(IntendedError::HATCH, -n, n);

    // GPSHandler handler = GPSHandler("brdc0600.26n");
    // SatNav problem = SatNav("20260301", SWARM::A, handler);
    // SatNav problem = SatNav("20260301", SWARM::B, handler);
    // SatNav problem = SatNav("20260301", SWARM::C, handler);

    // Plotter plotter(problem, 0, 86400);
    // plotter.plot_errors_pr(-10, 10);
    // plotter.plot_errors_norm(0, 10);

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

    double m = 5;
    plotter.plot_errors_norm(0, m);
    plotter.plot_errors_proj(-m, m);
    plotter.plot_true_norm();

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