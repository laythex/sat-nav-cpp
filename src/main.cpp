#include "GPSHandler.hpp"
#include "SatNav.hpp"
#include "Plotter.hpp"

int main() {
    // GPSHandler handler = GPSHandler("brdc2940.19n"); SatNav problem = SatNav("GNV1B_2019-10-21_C_04.txt", "GPS1B_2019-10-21_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc1410.20n"); SatNav problem = SatNav("GNV1B_2020-05-20_C_04.txt", "GPS1B_2020-05-20_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc0890.21n"); SatNav problem = SatNav("GNV1B_2021-03-30_C_04.txt", "GPS1B_2021-03-30_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc2610.22n"); SatNav problem = SatNav("GNV1B_2022-09-18_C_04.txt", "GPS1B_2022-09-18_C_04.txt", handler);
    GPSHandler handler = GPSHandler("brdc0820.23n"); SatNav problem = SatNav("GNV1B_2023-03-23_C_04.txt", "GPS1B_2023-03-23_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc0940.24n"); SatNav problem = SatNav("GNV1B_2024-04-03_C_04.txt", "GPS1B_2024-04-03_C_04.txt", handler);

    Plotter plotter(problem);
    plotter.plot_errors_norm(0, 30);
    plotter.plot_errors_pr();

    return 0;
}