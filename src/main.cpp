#include "GPSHandler.hpp"
#include "NavigationalProblem.hpp"
#include "LinearAlgebra.hpp"

int main() {
    unsigned ti = 732801600 + 3600;
    unsigned tf = 732801600 + 3600 + 5400;

    std::cout << std::setprecision(15);
    GPSHandler handler = GPSHandler("../data/brdc0820.23n");
    NavigationalProblem problem = NavigationalProblem("../data/GNV1B_2023-03-23_C_04.txt", 
                                                      "../data/GPS1B_2023-03-23_C_04.txt", 
                                                      handler);
    problem.solve();
    problem.errors_norm("../results/errors-norm.txt");
    problem.errors_prs("../results/errors-prs.txt");

    return 0;
}