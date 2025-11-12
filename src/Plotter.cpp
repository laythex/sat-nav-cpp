#include "Plotter.hpp"

Plotter::Plotter(SatNav& problem) : problem(problem) {
    problem.solve();
}

void Plotter::plot_errors_norm(double ymin, double ymax) {
    std::string filename = "errors-norm";

    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Модуль ошибки" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        unsigned time = ss.time;
        State ts = problem.get_true_state_at(time);

        double norm = abs(ss.position - ts.position);

        file << time << sep << norm << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void Plotter::plot_errors_prs() {

}

void Plotter::run_py_plotter(const std::string& arg) const {
    std::string command = "python3 ../scripts/plotter.py " + arg;
    system(command.c_str());
}