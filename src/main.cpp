#include "GPSHandler.hpp"

int main() {

    GPSHandler handler = GPSHandler("../data/brdc0820.23n");
    std::vector<double> state = handler.interp_state(13, 732801618 + 90 * 60);
    std::cout << std::setprecision(10);
    std::cout << state[0] << '\t' << state[1] << '\t' << state[2] << '\t' << state[3] << '\t' << std::endl;

    return 0;
}