#include "Plotter.hpp"

Plotter::Plotter(const SatNav& sn, unsigned ti, unsigned tf) : problem(sn), problem_copy(sn), ti(ti), tf(tf) {
    problem.solve('0', ti, tf);
}

void Plotter::plot_errors_norm(double ymin, double ymax) {
    std::string filename = "errors-norm";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Модуль ошибки" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        if (!ss.is_solved) continue;

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

void Plotter::plot_errors_by_type(char error_type, double ymin, double ymax) {
    problem_copy.solve('I', ti, tf);

    std::string filename = "errors-" + error_names[error_type];
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << error_titles[error_type] << '\t' << "Время, с" << '\t' << "Вклад, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        if (!ss.is_solved) continue;

        unsigned time = ss.time;
        SolutionState copy_ss = problem_copy.get_solution_state_at(time);
        if (!copy_ss.is_solved) continue;

        double error_norm = abs(ss.position - copy_ss.position);

        file << time << sep << error_norm << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void Plotter::plot_map_norm(const std::string& mode, double threshold) {
    problem_copy.solve('I', ti, tf);

    std::string filename = "map-norm";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Ошибка модуля решения" << '\t' << mode << std::endl;
    file << threshold << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        if (!ss.is_solved) continue;

        unsigned time = ss.time;
        State ts = problem.get_true_state_at(time);

        std::vector<double> geo_coords = ECEF_to_geographycal(ss.position);
        double error_norm = log(abs(ss.position - ts.position));

        file << geo_coords[0] << sep << geo_coords[1] << sep << error_norm << std::endl;
    }

    file.close();

    run_py_mapper(filename);
}

void Plotter::plot_map_iono(const std::string& mode, double threshold) {
    problem_copy.solve('I', ti, tf);

    std::string filename = "map-iono";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Ионосферная ошибка" << '\t' << mode << std::endl;
    file << threshold << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        if (!ss.is_solved) continue;

        unsigned time = ss.time;
        SolutionState copy_ss = problem_copy.get_solution_state_at(time);
        if (!copy_ss.is_solved) continue;

        std::vector<double> geo_coords = ECEF_to_geographycal(ss.position);
        double error_norm = abs(ss.position - copy_ss.position);

        file << geo_coords[0] << sep << geo_coords[1] << sep << error_norm << std::endl;
    }

    file.close();

    run_py_mapper(filename);
}

void Plotter::run_py_plotter(const std::string& arg) const {
    std::string command = "python3 ../scripts/plotter.py " + arg;
    system(command.c_str());
}

void Plotter::run_py_mapper(const std::string& arg) const {
    std::string command = "python3 ../scripts/mapper.py " + arg;
    system(command.c_str());
}

std::vector<double> Plotter::ECEF_to_geographycal(const std::vector<double>& position) {
    double latitude, longitude;

    double x = position[0];
    double y = position[1];
    double z = position[2];
    double r = sqrt(x * x + y * y + z * z);

    longitude = atan2(x, y) * 180.0 * M_1_PI;
    latitude = asin(z / r) * 180.0 * M_1_PI;

    return {longitude, latitude};  
}