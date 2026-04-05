#include "GPSHandler.hpp"
#include "SatNav.hpp"
#include "Plotter.hpp"
#include "PlotterRel.hpp"

#include "Propagator.hpp"


int main() {

    // Разобраться с дельтами
    // Починить SNR в 2005
    // PRN_ID в struct?

    // GPSHandler handler = GPSHandler("brdc3440.05n"); SatNav problem = SatNav("2005-12-10", GRACE_SATS::A, handler);
    // GPSHandler handler = GPSHandler("brdc3440.05n"); SatNav problem = SatNav("2005-12-10", GRACE_SATS::B, handler);
    // GPSHandler handler = GPSHandler("brdc2940.19n"); SatNav problem = SatNav("2019-10-21", GRACE_SATS::C, handler);
    // GPSHandler handler = GPSHandler("brdc1410.20n"); SatNav problem = SatNav("2020-05-20", GRACE_SATS::C, handler);
    // GPSHandler handler = GPSHandler("brdc0890.21n"); SatNav problem = SatNav("2021-03-30", GRACE_SATS::C, handler);
    // GPSHandler handler = GPSHandler("brdc2610.22n"); SatNav problem = SatNav("2022-09-18", GRACE_SATS::C, handler);
    // GPSHandler handler = GPSHandler("brdc0820.23n"); SatNav problem = SatNav("2023-03-23", GRACE_SATS::C, handler);
    // GPSHandler handler = GPSHandler("brdc0940.24n"); SatNav problem = SatNav("2024-04-03", GRACE_SATS::C, handler);

    // Plotter plotter(problem, 0, 43200);
    // plotter.plot_errors_pr(-10, 10);
    // plotter.plot_errors_norm(0, 10);


    GPSHandler handler = GPSHandler("brdc0600.26n");
    SatNav problem = SatNav("20260301", SWARM_SATS::A, handler);
    // SatNav problem = SatNav("20260301", SWARM_SATS::B, handler);
    // SatNav problem = SatNav("20260301", SWARM_SATS::C, handler);

    Plotter plotter(problem, 0, 86400);
    plotter.plot_errors_pr(-10, 10);
    plotter.plot_errors_norm(0, 10);

    // GPSHandler handler = GPSHandler("brdc3440.05n");
    // SatNav problem1 = SatNav("2005-12-10", GRACE_SATS::A, '2', handler);
    // SatNav problem2 = SatNav("2005-12-10", GRACE_SATS::B, '2', handler);

    // SatNavRel problem_rel(problem1, problem2);
    // PlotterRel plotter(problem_rel, 0, 0);

    // double m = 5;
    // plotter.plot_errors_norm(0, m);
    // plotter.plot_errors_proj(-m, m);
    // plotter.plot_true_norm(0, 0);

    return 0;
}