#include <fstream>

#include "DataParser.hpp"
#include "PropagatorRel.hpp"

int main() {
    Date date1 = {2023, 3, 23};
    std::vector<State> ts1_pas = DataParser::load_grace_fo_gnv_data(date1, GRACE::C);
    std::vector<State> ts1_act = DataParser::load_grace_fo_gnv_data(date1, GRACE::D);

    Date date2 = {2005, 12, 10};
    std::vector<State> ts2_pas = DataParser::load_grace_gnv_data(date2, GRACE::A);
    std::vector<State> ts2_act = DataParser::load_grace_gnv_data(date2, GRACE::B);

    std::vector<State> ts1_rel(ts1_pas.size());
    for (size_t i = 0; i < ts1_pas.size(); ++i) {
        ts1_rel[i] = {ts1_pas[i].time, ts1_act[i].position - ts1_pas[i].position,
                      ts1_act[i].velocity - ts1_pas[i].velocity};
    }

    std::vector<State> ts2_rel(ts2_pas.size());
    for (size_t i = 0; i < ts2_pas.size(); ++i) {
        ts2_rel[i] = {ts2_pas[i].time, ts2_act[i].position - ts2_pas[i].position,
                      ts2_act[i].velocity - ts2_pas[i].velocity};
    }

    size_t i0 = 10;
    size_t j0 = 2;

    State i1_rk4 = ts1_rel[i0];
    State i1_inc = {ts1_rel[i0].time, ts1_rel[i0].position, ts1_rel[i0].position - ts1_rel[0].position};

    State i2_rk4 = ts2_rel[j0];
    State i2_inc = {ts2_rel[j0].time, ts2_rel[j0].position, ts2_rel[j0].position - ts2_rel[0].position};

    PropagatorRel p1_rk4(i1_rk4);
    PropagatorRel p1_inc(i1_inc);
    PropagatorRel p2_rk4(i2_rk4);
    PropagatorRel p2_inc(i2_inc);

    char sep = ',';
    std::ofstream file;
    file.open("../results/propagators-rel.csv", std::fstream::out);

    file << "Ошибка интегрирования" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;
    file << 0.0 << '\t' << 0.0 << std::endl;
    file << "Метод Рунге-Кутты, групповой полет";
    file << '\t' << "В приращениях, групповой полет";
    file << '\t' << "Метод Рунге-Кутты, сближение и стыковка";
    file << '\t' << "В приращениях, сближение и стыковка";
    file << std::endl;

    file << ts1_rel[i0].time - ts1_rel[i0].time;
    file << sep << (ts1_rel[i0].position - i1_rk4.position).norm();
    file << sep << (ts1_rel[i0].position - i1_inc.position).norm();
    file << sep << (ts2_rel[j0].position - i2_rk4.position).norm();
    file << sep << (ts2_rel[j0].position - i2_inc.position).norm();
    file << std::endl;

    size_t i = i0;
    size_t j = j0;
    for (size_t t = 10; t < 70; t += 10) {
        State s1_rk4 = p1_rk4.step_rk4_rel(10.0, ts1_pas[i]);
        State s1_inc = p1_inc.step_inc_rel(10.0, ts1_pas[i]);
        State s2_rk4 = p2_rk4.step_rk4_rel(10.0, ts2_pas[j]);
        State s2_inc = p2_inc.step_inc_rel(10.0, ts2_pas[j]);

        i += 10;
        j += 2;

        file << ts1_rel[i].time - ts1_rel[i0].time;
        file << sep << (ts1_rel[i].position - s1_rk4.position).norm();
        file << sep << (ts1_rel[i].position - s1_inc.position).norm();
        file << sep << (ts2_rel[j].position - s2_rk4.position).norm();
        file << sep << (ts2_rel[j].position - s2_inc.position).norm();
        file << std::endl;
    }

    file.close();

    std::string command = "python3 ../scripts/plotter.py propagators-rel";
    system(command.c_str());

    return 0;
}