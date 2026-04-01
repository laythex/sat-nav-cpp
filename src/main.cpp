#include "GPSHandler.hpp"
#include "SatNav.hpp"
#include "Plotter.hpp"
#include "PlotterRel.hpp"

#include "Propagator.hpp"


int main() {
    // GPSHandler handler = GPSHandler("brdc3440.05n"); SatNav problem = SatNav("GNV1B_2005-12-10_A_02.dat", "GPS1B_2005-12-10_A_02.dat", handler);
    // GPSHandler handler = GPSHandler("brdc3440.05n"); SatNav problem = SatNav("GNV1B_2005-12-10_B_02.dat", "GPS1B_2005-12-10_B_02.dat", handler);
    // GPSHandler handler = GPSHandler("brdc2940.19n"); SatNav problem = SatNav("GNV1B_2019-10-21_C_04.txt", "GPS1B_2019-10-21_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc1410.20n"); SatNav problem = SatNav("GNV1B_2020-05-20_C_04.txt", "GPS1B_2020-05-20_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc0890.21n"); SatNav problem = SatNav("GNV1B_2021-03-30_C_04.txt", "GPS1B_2021-03-30_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc2610.22n"); SatNav problem = SatNav("GNV1B_2022-09-18_C_04.txt", "GPS1B_2022-09-18_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc0820.23n"); SatNav problem = SatNav("GNV1B_2023-03-23_C_04.txt", "GPS1B_2023-03-23_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc0940.24n"); SatNav problem = SatNav("GNV1B_2024-04-03_C_04.txt", "GPS1B_2024-04-03_C_04.txt", handler);

    // Plotter plotter(problem, 0, 0);
    // plotter.plot_errors_norm(0, 20);

    // GPSHandler handler = GPSHandler("brdc3450.05n"); // ИСПРАВИТЬ
    // SatNav problem1 = SatNav("2005-12-11", SatNav::GRACE_SATS::A, '2', handler);
    // SatNav problem2 = SatNav("2005-12-11", SatNav::GRACE_SATS::B, '2', handler);

    // GPSHandler handler = GPSHandler("brdc3450.05n");
    // SatNav problem1 = SatNav("20260301", SatNav::SWARM_SATS::A, handler);
    // SatNav problem2 = SatNav("20260301", SatNav::SWARM_SATS::B, handler);
    // SatNav problem3 = SatNav("20260301", SatNav::SWARM_SATS::C, handler);

    DataParser::load_swarm_nav_data("SW_OPER_GPSANAV_1B_20260301T000000_20260301T235959_0602.sp3");

    // SatNavRel problem_rel(problem1, problem2);
    // PlotterRel plotter(problem_rel, 0, 0);

    // double m = 5;
    // plotter.plot_errors_norm(0, m);
    // plotter.plot_errors_proj(-m, m);
    // plotter.plot_true_norm(0, 0);

    return 0;
}