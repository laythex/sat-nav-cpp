#include "GPSHandler.hpp"
#include "NavigationalProblem.hpp"
#include "Matrix.hpp"

int main() {
    
    // GPSHandler handler = GPSHandler("../data/brdc0820.23n");
    // NavigationalProblem problem = NavigationalProblem("../data/GNV1B_2023-03-23_C_04.txt", 
    //                                                   "../data/GPS1B_2023-03-23_C_04.txt", 
    //                                                   handler);
    // problem.solve(732801610, 732801695);

    Matrix mat({{1, 2, 3}, {3, 1, 2}, {2, 3, 1}});

    std::cout << mat << std::endl;
    std::cout << mat.transpose() << std::endl;

    return 0;
}