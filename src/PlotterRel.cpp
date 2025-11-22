#include "PlotterRel.hpp"

PlotterRel::PlotterRel(const SatNavRel& sn, unsigned ti, unsigned tf) : problem(sn), ti(ti), tf(tf) {
    problem.solve_relative('0', ti, tf);
}

void PlotterRel::plot_errors_norm(double ymin, double ymax) {
    std::string filename = "rel-norm";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Модуль расстояния" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        if (!ss.is_solved) continue;

        unsigned time = ss.time;

        double ss_norm = abs(ss.position);

        // переписать через тру стейт проблемы
        double ts_norm = abs(problem.pas.get_true_state_at(time).position - problem.act.get_true_state_at(time).position);

        file << time << sep << ss_norm << sep << ts_norm << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void PlotterRel::run_py_plotter(const std::string& arg) const {
    std::string command = "python3 ../scripts/plotter.py " + arg;
    system(command.c_str());
}

std::vector<double> PlotterRel::ECEF_to_geographycal(const std::vector<double>& position) {
    double latitude, longitude;

    double x = position[0];
    double y = position[1];
    double z = position[2];
    double r = sqrt(x * x + y * y + z * z);

    longitude = atan2(x, y) * 180.0 * M_1_PI;
    latitude = asin(z / r) * 180.0 * M_1_PI;

    return {longitude, latitude};  
}