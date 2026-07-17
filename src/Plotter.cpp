#include "Plotter.hpp"

Plotter::Plotter(const SatNav& sn) : problem(sn), problem_copy(sn) { problem.solve(IntendedError::NONE); }

void Plotter::plot_errors_norm(double ymin, double ymax) {
    std::string filename = "errors-norm-" + to_string(problem.get_date());
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Модуль ошибки" << '\t' << "Время, ч" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << std::endl;

    unsigned t0 = problem.get_solution_states().front().time;
    for (const auto& ss : problem.get_solution_states()) {
        if (ss.status != SolutionStatus::OK) continue;

        double time = (ss.time - t0) / 3600.0;

        auto ts_it = problem.get_true_state_iterator(ss.time);
        if (ts_it == problem.get_true_states().end()) continue;
        State ts = *ts_it;

        double error_norm = (ss.position - ts.position).norm();

        file << time << sep << error_norm << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void Plotter::plot_errors_pr(double ymin, double ymax) {
    std::string filename = "errors-pr-" + to_string(problem.get_date());
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Модуль ошибки псевдодальности" << '\t' << "Время, ч" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << std::endl;

    double dt_ASN_c_avg = 0.0;
    unsigned ss_count = 0;

    for (const auto& ss : problem.get_solution_states()) {
        if (ss.status != SolutionStatus::OK) continue;

        dt_ASN_c_avg += ss.dt_ASN_c;
        ++ss_count;
    }

    dt_ASN_c_avg /= ss_count;

    unsigned t0 = problem.get_refined_measurements_groupped().front().time;
    for (const auto& ref_mg : problem.get_refined_measurements_groupped()) {
        double time = (ref_mg.time - t0) / 3600.0;

        auto ts_it = problem.get_true_state_iterator(ref_mg.time);
        if (ts_it == problem.get_true_states().end()) continue;
        State ts = *ts_it;

        auto ss_it = problem.get_solution_state_iterator(ref_mg.time);
        if (ss_it == problem.get_solution_states().end()) continue;
        SolutionState ss = *ss_it;

        file << time << sep;

        for (unsigned prn_id = 1; prn_id <= 32; ++prn_id) {
            unsigned prn_index = prn_id - 1;
            RefinedMeasurement ref_m = ref_mg.measurements[prn_index];

            if (ref_m.status == MeasurementStatus::OK && ss.status == SolutionStatus::OK) {
                double error_pr =
                    (ss.gps_states[prn_index].position - ts.position).norm() - ref_m.pr_refined - dt_ASN_c_avg;
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

void Plotter::plot_errors_norm_by_type(IntendedError error, double ymin, double ymax) {
    SatNav sn = problem_copy;
    sn.solve(error);

    std::string filename = "errors-norm-" + to_filename(error) + "-" + to_string(problem.get_date());
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Исследуемая ошибка: " + to_title(error) << '\t' << "Время, ч" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << std::endl;

    unsigned t0 = problem.get_solution_states().front().time;
    for (const auto& ss : problem.get_solution_states()) {
        if (ss.status != SolutionStatus::OK) continue;

        unsigned time = ss.time;

        auto ss_copy_it = sn.get_solution_state_iterator(time);
        if (ss_copy_it == sn.get_solution_states().end()) continue;
        SolutionState ss_copy = *ss_copy_it;

        if (ss_copy.status != SolutionStatus::OK) continue;

        State ts = *problem.get_true_state_iterator(time);

        // double error_norm = (ss_copy.position - ss.position).norm();

        file << (time - t0) / 3600.0 << sep << (ss.position - ts.position).norm() << sep
             << (ss_copy.position - ts.position).norm() << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void Plotter::plot_errors_proj_by_type(IntendedError error, double ymin, double ymax) {
    SatNav sn = problem_copy;
    sn.solve(error);

    std::string filename = "errors-proj-" + to_filename(error) + "-" + to_string(problem.get_date());
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Вклад в решение: " + to_title(error) << '\t' << "Время, ч" << '\t' << "Вклад, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << 'x' << '\t' << 'y' << '\t' << 'z' << std::endl;

    unsigned t0 = problem.get_solution_states().front().time;
    for (const auto& ss : problem.get_solution_states()) {
        if (ss.status != SolutionStatus::OK) continue;

        double time = (ss.time - t0) / 3600.0;

        auto ss_copy_it = sn.get_solution_state_iterator(ss.time);
        if (ss_copy_it == sn.get_solution_states().end()) continue;
        SolutionState ss_copy = *ss_copy_it;

        if (ss_copy.status != SolutionStatus::OK) continue;

        Eigen::Vector3d error = ss_copy.position - ss.position;

        file << time << sep << error[0] << sep << error[1] << sep << error[2] << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void Plotter::plot_errors_pr_by_type(IntendedError error, double ymin, double ymax) {
    SatNav sn = problem_copy;
    sn.solve(error);

    std::string filename = "errors-pr-" + to_filename(error) + "-" + to_string(problem.get_date());
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Ошибка псевдодальности без: " + to_title(error) << '\t' << "Время, ч" << '\t' << "Ошибка, м" << std::endl;
    file << ymin << '\t' << ymax << std::endl;
    file << std::endl;

    double dt_ASN_c_avg = 0.0;
    unsigned ss_count = 0;

    for (const auto& ss_copy : sn.get_solution_states()) {
        if (ss_copy.status != SolutionStatus::OK) continue;

        dt_ASN_c_avg += ss_copy.dt_ASN_c;
        ++ss_count;
    }

    dt_ASN_c_avg /= ss_count;

    unsigned t0 = problem.get_refined_measurements_groupped().front().time;
    for (const auto& ref_mg : sn.get_refined_measurements_groupped()) {
        double time = (ref_mg.time - t0) / 3600.0;

        auto ts_copy_it = sn.get_true_state_iterator(ref_mg.time);
        if (ts_copy_it == sn.get_true_states().end()) continue;
        State ts = *ts_copy_it;

        auto ss_copy_it = sn.get_solution_state_iterator(ref_mg.time);
        if (ss_copy_it == sn.get_solution_states().end()) continue;
        SolutionState ss_copy = *ss_copy_it;

        file << time << sep;

        for (unsigned prn_id = 1; prn_id <= 32; ++prn_id) {
            unsigned prn_index = prn_id - 1;
            RefinedMeasurement ref_m = ref_mg.measurements[prn_index];

            if (ref_m.status == MeasurementStatus::OK && ss_copy.status == SolutionStatus::OK) {
                double error_pr =
                    (ss_copy.gps_states[prn_index].position - ts.position).norm() - ref_m.pr_refined - dt_ASN_c_avg;
                file << error_pr;
            }

            file << sep;
        }

        file << std::endl;
    }

    file.close();

    run_py_plotter(filename);
}

void Plotter::plot_map_iono(const std::string& mode, double threshold) {
    SatNav sn = problem_copy;
    sn.solve(IntendedError::IONOSPHERIC);

    std::string filename = "map-iono";
    std::ofstream file;
    file.open("../results/" + filename + ".csv", std::fstream::out);

    file << "Ионосферная ошибка" << '\t' << mode << std::endl;
    file << threshold << std::endl;

    for (const auto& ss : problem.get_solution_states()) {
        if (ss.status != SolutionStatus::OK) continue;

        unsigned time = ss.time;

        auto ss_copy_it = sn.get_solution_state_iterator(time);
        if (ss_copy_it == sn.get_solution_states().end()) continue;
        SolutionState ss_copy = *ss_copy_it;

        if (ss_copy.status != SolutionStatus::OK) continue;

        Eigen::Vector3d geo_coords = ECEF_to_geographycal(ss.position);
        double error_norm = (ss.position - ss_copy.position).norm();

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
        file << acc_m.time << sep << acc_m.linear_acceleration[0] << sep << acc_m.linear_acceleration[1] << sep
             << acc_m.linear_acceleration[2] << std::endl;
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
        file << acc_m.time << sep << acc_m.angular_acceleration[0] << sep << acc_m.angular_acceleration[1] << sep
             << acc_m.angular_acceleration[2] << std::endl;
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

Eigen::Vector3d Plotter::ECEF_to_geographycal(const Eigen::Vector3d& position) {
    double latitude, longitude;

    double x = position[0];
    double y = position[1];
    double z = position[2];
    double r = sqrt(x * x + y * y + z * z);

    longitude = atan2(x, y) * 180.0 * M_1_PI;
    latitude = asin(z / r) * 180.0 * M_1_PI;

    return {longitude, latitude, 0.0};
}