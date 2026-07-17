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
    Config::CN0_min_threshold = 25.0;
    Config::CN0_max_threshold = 60.0;
    Config::fadein_threshold = 60;
    Config::hatch_constant = 1 / 10.0;
    Config::solution_tolerance = 0.1;
    Config::residual_threshold = 1.5;
    Config::residual_rel_threshold = 10.0;

    // Date date = {2005, 12, 10};
    // Date date = {2005, 12, 10};
    // Date date = {2019, 10, 21};
    // Date date = {2020, 5, 20};
    // Date date = {2021, 3, 30};
    // Date date = {2022, 9, 18};
    Date date = {2023, 3, 23};
    // Date date = {2024, 4, 3};
    // Date date = {2024, 10, 1};

    GRACE sat = GRACE::C;

    GPSHandler handler = GPSHandler(date);
    SatNav problem = SatNav(date, sat, handler);

    double m = 10;
    Plotter plotter(problem);
    plotter.plot_errors_norm(0, m);

    // plotter.plot_errors_norm_by_type(IntendedError::RELATIVISTIC, 0.0, 30.0);
    // plotter.plot_errors_norm_by_type(IntendedError::IONOSPHERIC, 0.0, 30.0);
    // plotter.plot_errors_norm_by_type(IntendedError::SAGNAC, 0.0, 50.0);
    // plotter.plot_errors_norm_by_type(IntendedError::NOISY, 0.0, 10.0);
    // plotter.plot_errors_norm_by_type(IntendedError::HATCH, 0.0, 10.0);

    // GPSHandler handler = GPSHandler("brdc0600.26n");
    // SatNav problem = SatNav("20260301", SWARM::A, handler);
    // SatNav problem = SatNav("20260301", SWARM::B, handler);
    // SatNav problem = SatNav("20260301", SWARM::C, handler);

    // Plotter plotter(problem, 0, 86400);
    // plotter.plot_errors_pr(-10, 10);
    // plotter.plot_errors_norm(0, 10);

    return 0;
}