#include "GPSHandler.hpp"
#include "SatNav.hpp"
#include "LinAlg.hpp"

int main() {

    // GPSHandler handler = GPSHandler("brdc0820.23n");
    // SatNav problem = SatNav("GNV1B_2023-03-23_C_04.txt", "GPS1B_2023-03-23_C_04.txt", handler);
    GPSHandler handler = GPSHandler("brdc0940.24n");
    SatNav problem = SatNav("GNV1B_2024-04-03_C_04.txt", "GPS1B_2024-04-03_C_04.txt", handler);

    problem.solve();
    problem.out_error_norm();
    problem.out_error_prs();
    // problem.out_number_of_sats();
    // problem.out_snr_over_sol();

    // for (unsigned i = 0; i < 5; i++) {
    //     problem.out_error_by_type(i);
    // }

    return 0;
}