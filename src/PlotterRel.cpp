#include "PlotterRel.hpp"

PlotterRel::PlotterRel(const SatNavRel& sn) : problem(sn), problem_copy(sn) { problem.solve(IntendedError::NONE); }

void PlotterRel::plot_errors_norm(double ymin, double ymax) {
    std::string filename = "rel-errors-norm-" + to_string(problem.pas.get_date());
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Модуль ошибки" << '\t' << "Время, ч" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << "Измерения" << '\t' << "Модель" << std::endl;

    unsigned t0 = problem.get_solution_states().front().time;
    for (const auto& ss : problem.get_solution_states()) {
        if (ss.status != SolutionStatusRel::OK && ss.status != SolutionStatusRel::MODEL_STEP) continue;

        double time = (ss.time - t0) / 3600.0;

        auto ts_it = problem.get_true_state_iterator(ss.time);
        if (ts_it == problem.get_true_states().end()) continue;
        State ts = *ts_it;

        double error_norm = (ts.position - ss.position).norm();

        if (ss.status == SolutionStatusRel::OK) {
            file << time << sep << error_norm << sep << -1.0 << std::endl;
        } else if (ss.status == SolutionStatusRel::MODEL_STEP) {
            file << time << sep << -1.0 << sep << error_norm << std::endl;
        }
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
        if (ss.status != SolutionStatusRel::OK) continue;

        unsigned time = ss.time;

        auto ts_it = problem.get_true_state_iterator(time);
        if (ts_it == problem.get_true_states().end()) continue;
        State ts = *ts_it;

        Eigen::Vector3d error = ts.position - ss.position;

        file << time << sep << error[0] << sep << error[1] << sep << error[2] << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void PlotterRel::plot_true_norm(double ymin, double ymax) {
    std::string filename = "true-norm-" + to_string(problem.pas.get_date());
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Расстояние между аппаратами" << '\t' << "Время, ч" << '\t' << "Расстояние, км" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << std::endl;

    unsigned t0 = problem.get_true_states().front().time;
    for (const auto& ts : problem.get_true_states()) {
        double time = (ts.time - t0) / 3600.0;
        double true_norm = ts.position.norm();

        file << time << sep << true_norm * 1.0e-3 << std::endl;
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

void PlotterRel::plot_errors_norm_by_type(IntendedError error, double ymin, double ymax) {
    SatNavRel sn = problem_copy;
    sn.solve(error);

    std::string filename = "rel-errors-norm-" + to_filename(error) + "-" + to_string(problem.pas.get_date());
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Исследуемая ошибка: " + to_title(error) << '\t' << "Время, ч" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << std::endl;

    unsigned t0 = problem.get_solution_states().front().time;
    for (const auto& ss : problem.get_solution_states()) {
        if (ss.status != SolutionStatusRel::OK) continue;

        unsigned time = ss.time;

        auto ss_copy_it = sn.get_solution_state_iterator(time);
        if (ss_copy_it == sn.get_solution_states().end()) continue;
        SolutionStateRel ss_copy = *ss_copy_it;

        if (ss_copy.status != SolutionStatusRel::OK) continue;

        State ts = *problem.get_true_state_iterator(time);

        // double error_norm = (ss_copy.position - ss.position).norm();

        file << (time - t0) / 3600.0 << sep << (ss.position - ts.position).norm() << sep
             << (ss_copy.position - ts.position).norm() << std::endl;
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