#include "GPSHandler.hpp"
#include "NavigationalProblem.hpp"

// unordered_map?

int main() {
    
    GPSHandler handler = GPSHandler("../data/brdc0820.23n");
    NavigationalProblem problem = NavigationalProblem("../data/GNV1B_2023-03-23_C_04.txt", "../data/GPS1B_2023-03-23_C_04.txt");
    problem.solve(732801610, 732801695);
    // std::cout << handler.get_state(13, 732456018 + 345600)[1] << std::endl;

    return 0;
}