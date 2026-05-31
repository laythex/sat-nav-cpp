#include "GPSHandler.hpp"
#include "SatNav.hpp"
#include "Plotter.hpp"
#include "PlotterRel.hpp"

int main() {
    /*
    TODO:
    все string через format
    PRN_ID в struct?
    Raw Measurement хранит C1C, L1C и т.д.
    Error type и failure type в enum
    */

    // Date date = {2005, 12, 10}; GRACE sat = GRACE::A;
    // Date date = {2005, 12, 10}; GRACE sat = GRACE::B;
    // Date date = {2019, 10, 21}; GRACE sat = GRACE::C;
    // Date date = {2020, 5, 20}; GRACE sat = GRACE::C;
    // Date date = {2021, 3, 30}; GRACE sat = GRACE::C;
    // Date date = {2022, 9, 18}; GRACE sat = GRACE::C;
    // Date date = {2023, 3, 23}; GRACE sat = GRACE::C;
    Date date = {2024, 4, 3}; GRACE sat = GRACE::C;

    GPSHandler handler = GPSHandler(date); 
    SatNav problem = SatNav(date, sat, handler);

    double m = 10;
    Plotter plotter(problem, 0, 0);
    plotter.plot_errors_pr(-m, m);
    plotter.plot_errors_norm(0, m);

    // GPSHandler handler = GPSHandler("brdc0600.26n");
    // SatNav problem = SatNav("20260301", SWARM::A, handler);
    // SatNav problem = SatNav("20260301", SWARM::B, handler);
    // SatNav problem = SatNav("20260301", SWARM::C, handler);

    // Plotter plotter(problem, 0, 86400);
    // plotter.plot_errors_pr(-10, 10);
    // plotter.plot_errors_norm(0, 10);

    // Date date = {2005, 12, 9};
    // GPSHandler handler = GPSHandler(date);
    // SatNav problem1 = SatNav(date, GRACE::A, handler);
    // SatNav problem2 = SatNav(date, GRACE::B, handler);

    // SatNavRel problem_rel(problem1, problem2);
    // PlotterRel plotter(problem_rel);

    // double m = 5;
    // plotter.plot_errors_norm(0, m);
    // plotter.plot_errors_proj(-m, m);
    // plotter.plot_true_norm();

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