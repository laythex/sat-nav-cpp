#include "GPSHandler.hpp"
#include "NavigationalProblem.hpp"
#include "LinearAlgebra.hpp"

int main() {

    GPSHandler handler = GPSHandler("brdc0820.23n");
    NavigationalProblem problem = NavigationalProblem("GNV1B_2023-03-23_C_04.txt", "GPS1B_2023-03-23_C_04.txt", handler);

    problem.solve();
    problem.out_errors_norm();
    problem.out_errors_prs();
    // problem.out_errors_rel();
    problem.out_number_of_sats();

    return 0;
}