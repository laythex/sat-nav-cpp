#include "PlotterRel.hpp"

PlotterRel::PlotterRel(const SatNavRel& sn, double ti, double tf) : problem(sn), ti(ti), tf(tf) {
    problem.solve_relative('0', ti, tf);
}

void PlotterRel::plot_errors_norm(double ymin, double ymax) {
    std::string filename = "rel-errors-norm";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Модуль ошибки" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        if (!ss.is_solved) continue;

        unsigned time = ss.time;

        auto ts_it = problem.get_true_state_iterator(time);
        if (ts_it == problem.get_true_states().end()) continue;
        State ts = *ts_it;

        double error_norm = abs(ts.position - ss.position);
        
        file << time << sep << error_norm << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void PlotterRel::plot_errors_proj(double ymin, double ymax) {
    std::string filename = "rel-errors-proj";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Проекция ошибки" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << 'x' << '\t' << 'y' << '\t' << 'z' << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        if (!ss.is_solved) continue;

        unsigned time = ss.time;

        auto ts_it = problem.get_true_state_iterator(time);
        if (ts_it == problem.get_true_states().end()) continue;
        State ts = *ts_it;

        std::vector<double> error = ts.position - ss.position;
        
        file << time << sep << error[0] << sep << error[1] << sep << error[2] << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void PlotterRel::plot_true_norm(double ymin, double ymax) {
    std::string filename = "true-norm";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Расстояние между аппаратами" << '\t' << "Время, с" << '\t' << "Расстояние, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << std::endl;

    for (const auto& ts : problem.get_true_states()) {
        unsigned time = ts.time;
        double true_norm = abs(ts.position);
        
        file << time << sep << true_norm << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void PlotterRel::plot_true_proj(double ymin, double ymax) { // Сделать легенду
    std::string filename = "true-proj";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Расстояние между аппаратами вдоль осей" << '\t' << "Время, с" << '\t' << "Расстояние, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << 'x' << '\t' << 'y' << '\t' << 'z' << std::endl;
    
    for (const auto& ts : problem.get_true_states()) {
        unsigned time = ts.time;
        
        file << time << sep << ts.position[0] << sep << ts.position[1] << sep << ts.position[2] << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void PlotterRel::run_py_plotter(const std::string& arg) const {
    std::string command = "python3 ../scripts/plotter.py " + arg;
    system(command.c_str());
}

std::vector<double> PlotterRel::ECEF_to_geographical(const std::vector<double>& position) {
    double latitude, longitude;

    double x = position[0];
    double y = position[1];
    double z = position[2];
    double r = sqrt(x * x + y * y + z * z);

    longitude = atan2(x, y) * 180.0 * M_1_PI;
    latitude = asin(z / r) * 180.0 * M_1_PI;

    return {longitude, latitude};  
}