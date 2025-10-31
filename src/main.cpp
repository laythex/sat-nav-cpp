#include "GPSHandler.hpp"
#include "NavigationalProblem.hpp"
#include "LinearAlgebra.hpp"

int main() {

    GPSHandler handler = GPSHandler("brdc0820.23n");
    NavigationalProblem problem = NavigationalProblem("GNV1B_2023-03-23_C_04.txt", "GPS1B_2023-03-23_C_04.txt", handler);

    problem.solve();
    problem.errors_norm();
    problem.errors_prs();
    problem.errors_rel();

    return 0;
}