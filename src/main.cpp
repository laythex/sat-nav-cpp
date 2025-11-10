#include "GPSHandler.hpp"
#include "SatNav.hpp"
#include "LinAlg.hpp"

int main() {
    // GPSHandler handler = GPSHandler("brdc2940.19n"); SatNav problem = SatNav("GNV1B_2019-10-21_C_04.txt", "GPS1B_2019-10-21_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc1410.20n"); SatNav problem = SatNav("GNV1B_2020-05-20_C_04.txt", "GPS1B_2020-05-20_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc0890.21n"); SatNav problem = SatNav("GNV1B_2021-03-30_C_04.txt", "GPS1B_2021-03-30_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc2610.22n"); SatNav problem = SatNav("GNV1B_2022-09-18_C_04.txt", "GPS1B_2022-09-18_C_04.txt", handler);
    GPSHandler handler = GPSHandler("brdc0820.23n"); SatNav problem = SatNav("GNV1B_2023-03-23_C_04.txt", "GPS1B_2023-03-23_C_04.txt", handler);
    // GPSHandler handler = GPSHandler("brdc0940.24n"); SatNav problem = SatNav("GNV1B_2024-04-03_C_04.txt", "GPS1B_2024-04-03_C_04.txt", handler);

    problem.solve();
    problem.out_error_norm();
    problem.out_error_prs();
    // problem.out_is_solved();
    // problem.out_number_of_sats();
    // problem.out_error_by_type(5);

    // for (unsigned i = 0; i < 5; i++) {
    //     problem.out_error_by_type(i);
    // }

    // std::map<unsigned, std::vector<double>> pos; 
    // std::map<unsigned, std::vector<double>> vel;
    // DataParser::load_gnv_data("GNV1B_2005-12-10_A_02.dat", pos, vel);

    return 0;
}