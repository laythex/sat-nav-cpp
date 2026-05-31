#include "Plotter.hpp"

Plotter::Plotter(const SatNav& sn, unsigned time_initial, unsigned time_final) : problem(sn), problem_copy(sn), ti(time_initial), tf(time_final) {
    problem.solve(ti, tf, IntendedError::NONE);
}

void Plotter::plot_errors_norm(double ymin, double ymax) {
    std::string filename = "errors-norm";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Модуль ошибки" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        if (ss.status != SolutionStatus::OK) continue;

        unsigned time = ss.time;

        auto ts_it = problem.get_true_state_iterator(time);
        if (ts_it == problem.get_true_states().end()) continue;
        State ts = *ts_it;

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
    file << std::endl;

    for (const auto& ref_mg : problem.get_refined_measurements_groupped()) {
        unsigned time = ref_mg.time;

        auto ts_it = problem.get_true_state_iterator(time);
        if (ts_it == problem.get_true_states().end()) continue;
        State ts = *ts_it;

        file << time << sep;

        for (unsigned prn_id = 1; prn_id <= 32; ++prn_id) {
            unsigned prn_index = prn_id - 1;
            RefinedMeasurement ref_m = ref_mg.refined_measurements[prn_index];

            if (ref_m.status == MeasurementStatus::OK) {
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

void Plotter::plot_gdop(double ymin, double ymax) {
std::string filename = "gdop";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "GDOP" << '\t' << "Время, с" << '\t' << "GDOP" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        if (ss.status != SolutionStatus::OK) continue;

        unsigned time = ss.time;

        double gdop = ss.GDOP;

        file << time << sep << gdop << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void Plotter::plot_errors_by_type(IntendedError error, double ymin, double ymax) {
    problem_copy.solve(ti, tf, error);

    std::string filename = "errors-" + to_filename(error);
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Исследуемая ошибка: " + to_title(error) << '\t' << "Время, с" << '\t' << "Вклад, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        if (ss.status != SolutionStatus::OK) continue;

        unsigned time = ss.time;

        auto copy_ss_it = problem.get_solution_state_iterator(time);
        if (copy_ss_it == problem.get_solution_states().end()) continue;
        SolutionState copy_ss = *copy_ss_it;

        if (copy_ss.status != SolutionStatus::OK) continue;

        double error_norm = abs(ss.position - copy_ss.position);

        file << time << sep << error_norm << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void Plotter::plot_map_norm(const std::string& mode, double threshold) {
    problem_copy.solve(ti, tf, IntendedError::IONOSPHERIC);

    std::string filename = "map-norm";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Ошибка модуля решения" << '\t' << mode << std::endl;
    file << threshold << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        if (ss.status != SolutionStatus::OK) continue;

        unsigned time = ss.time;

        auto ts_it = problem.get_true_state_iterator(time);
        if (ts_it == problem.get_true_states().end()) continue;
        State ts = *ts_it;

        std::vector<double> geo_coords = ECEF_to_geographycal(ss.position);
        double error_norm = log(abs(ss.position - ts.position));

        file << geo_coords[0] << sep << geo_coords[1] << sep << error_norm << std::endl;
    }

    file.close();

    run_py_mapper(filename);
}

void Plotter::plot_map_iono(const std::string& mode, double threshold) {
    problem_copy.solve(ti, tf, IntendedError::IONOSPHERIC);

    std::string filename = "map-iono";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Ионосферная ошибка" << '\t' << mode << std::endl;
    file << threshold << std::endl;
    
    for (const auto& ss : problem.get_solution_states()) {
        if (ss.status != SolutionStatus::OK) continue;

        unsigned time = ss.time;
        
        auto copy_ss_it = problem.get_solution_state_iterator(time);
        if (copy_ss_it == problem.get_solution_states().end()) continue;
        SolutionState copy_ss = *copy_ss_it;

        if (copy_ss.status != SolutionStatus::OK) continue;

        std::vector<double> geo_coords = ECEF_to_geographycal(ss.position);
        double error_norm = abs(ss.position - copy_ss.position);

        file << geo_coords[0] << sep << geo_coords[1] << sep << error_norm << std::endl;
    }

    file.close();

    run_py_mapper(filename);
}

void Plotter::plot_acc_lin_proj(double ymin, double ymax) {
    std::string filename = "acc-lin-proj";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Линейное ускорение вдоль осей" << '\t' << "Время, с" << '\t' << "Ускорение, м/c^2" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << 'x' << '\t' << 'y' << '\t' << 'z' << std::endl;

    for (const auto& acc_m : problem.get_acceleration_measurements()) {
        file << acc_m.time << sep << acc_m.linear_acceleration[0] << sep << 
                                     acc_m.linear_acceleration[1] << sep << 
                                     acc_m.linear_acceleration[2] << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void Plotter::plot_acc_ang_proj(double ymin, double ymax) {
    std::string filename = "acc-ang-proj";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Угловое ускорение вокруг осей" << '\t' << "Время, с" << '\t' << "Ускорение, рад/c^2" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << 'x' << '\t' << 'y' << '\t' << 'z' << std::endl;

    for (const auto& acc_m : problem.get_acceleration_measurements()) {
        file << acc_m.time << sep << acc_m.angular_acceleration[0] << sep << 
                                     acc_m.angular_acceleration[1] << sep << 
                                     acc_m.angular_acceleration[2] << std::endl;
    }

    file.close();

    run_py_plotter(filename);
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