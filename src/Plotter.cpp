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

        double error_norm = abs(ss.position - ts.position);

        file << time << sep << error_norm << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void Plotter::plot_errors_pr(double ymin, double ymax) {
    std::string filename = "errors-pr";

    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Модуль ошибки псевдодальности" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;

    for (const auto& ref_mg : problem.get_refined_measurements_groupped()) {
        unsigned time = ref_mg.time;
        State ts = problem.get_true_state_at(time);

        file << time << sep;

        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            unsigned prn_index = prn_id - 1;
            RefinedMeasurement ref_m = ref_mg.refined_measurements[prn_index];

            if (ref_m.is_present) {
                double error_pr = abs(ref_m.gps_position - ts.position) - ref_m.pseudorange;
                file << error_pr;
            }

            file << sep;
        }

        file << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void plot_errors_by_type(char error_type, double ymin = 0, double ymax = 0) {

}

void plot_map_iono() {
    
}

void Plotter::run_py_plotter(const std::string& arg) const {
    std::string command = "python3 ../scripts/plotter.py " + arg;
    system(command.c_str());
}