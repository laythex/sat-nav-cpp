#include <fstream>

#include "DataParser.hpp"
#include "Propagator.hpp"

int main() {
    Date date = {2023, 3, 23};
    std::vector<State> ts = DataParser::load_grace_fo_gnv_data(date, GRACE::C);

    size_t i0 = 10;

    State i_rk4 = ts[i0];
    State i_inc = {ts[i0].time, ts[i0].position, ts[i0].position - ts[0].position};

    Propagator p_rk4(i_rk4);
    Propagator p_inc(i_inc);
    Propagator p_rk4_j(i_rk4, false);
    Propagator p_inc_j(i_inc, false);

    char sep = ',';
    std::ofstream file;
    file.open("../results/propagators.csv", std::fstream::out);

    file << "Ошибка интегрирования" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;
    file << 0.0 << '\t' << 0.0 << std::endl;
    file << "Метод РК4";
    file << '\t' << "«В приращениях»";
    file << '\t' << "Метод РК4, без J2";
    file << '\t' << "«В приращениях», без J2";
    file << std::endl;

    file << ts[i0].time - ts[i0].time;
    file << sep << (ts[i0].position - i_rk4.position).norm();
    file << sep << (ts[i0].position - i_inc.position).norm();
    file << sep << (ts[i0].position - i_rk4.position).norm();
    file << sep << (ts[i0].position - i_inc.position).norm();
    file << std::endl;

    for (size_t i = i0; i < i0 + 60; i += i0) {
        State s_rk4 = p_rk4.step_rk4(10.0);
        State s_inc = p_inc.step_inc(10.0);
        State s_rk4_j = p_rk4_j.step_rk4(10.0);
        State s_inc_j = p_inc_j.step_inc(10.0);

        file << ts[i + i0].time - ts[i0].time;
        file << sep << (ts[i + i0].position - s_rk4.position).norm();
        file << sep << (ts[i + i0].position - s_inc.position).norm();
        file << sep << (ts[i + i0].position - s_rk4_j.position).norm();
        file << sep << (ts[i + i0].position - s_inc_j.position).norm();
        file << std::endl;
    }

    file.close();

    std::string command = "python3 ../scripts/plotter.py propagators";
    system(command.c_str());

    return 0;
}