#include "SatNavRel.hpp"
#include "SatNav.hpp"

SatNavRel::SatNavRel(SatNav& passive, SatNav& active) : pas(passive), act(active) {
    // Решаем каждую из задач отдельно
    pas.solve();
    act.solve();

    // Заполняем вектор истинных решений
    for (const State& ts_pas : pas.get_true_states()) {

        auto ts_act_it = act.get_true_state_iterator(ts_pas.time);
        if (ts_act_it == act.get_true_states().end()) continue;
        State ts_act = *ts_act_it;

        true_states.push_back({ts_pas.time, ts_act.position - ts_pas.position, ts_act.velocity - ts_pas.velocity});
    }
}

SatNavRel::SatNavRel(const SatNavRel& sn) : pas(sn.pas), act(sn.act),
                                            true_states(sn.true_states) {}

void SatNavRel::solve_relative(char et, int ti, int tf) {

    logger.log("Starting to solve...");

    error_type = et;

    for (const auto& raw_mg_pas : pas.raw_measurements_groupped) {

        if (ti > 0 || tf > 0) {
            if (raw_mg_pas.time < ti || raw_mg_pas.time >= tf) continue;
        }

        auto raw_mg_act_it = act.get_raw_measurement_groupped_iterator(raw_mg_pas.time);
        if (raw_mg_act_it == act.raw_measurements_groupped.end()) continue;
        RawMeasurementGroupped raw_mg_act = *raw_mg_act_it;

        if (raw_mg_act.time == -1) continue; // ???

        logger.log("Current time = " + std::to_string(raw_mg_pas.time));

        // Обрабатываем сырые измерения с каждого доступного спутника
        std::vector<RefinedMeasurement> ref_ms(32);
        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            unsigned prn_index = prn_id - 1;

            // Проверяем, что с измерениями все ок
            if (!check_raw(raw_mg_pas.raw_measurements[prn_index]) || !check_raw(raw_mg_act.raw_measurements[prn_index])) continue;

            logger.log("Refining " + std::to_string(prn_id));

            // Непосредственно их обрабатываем
            RefinedMeasurement ref_m = refine_raw(raw_mg_pas.raw_measurements[prn_index], raw_mg_act.raw_measurements[prn_index]);
            ref_ms[prn_index] = ref_m;
        }

        // Группируем обработанные измерения со всех спутников
        RefinedMeasurementGroupped ref_mg = {raw_mg_pas.time, ref_ms};
        refined_measurements_groupped.push_back(ref_mg);

        // Начинаем динамическую фильтрацию
        double lambda_x = T_x / (T_x + 1);
        double lambda_v = T_v / (T_v + 1);
        lambda.insert(identity(3) * lambda_x, 0, 0);
        lambda.insert(identity(3) * lambda_v, 3, 3);

        // Ищем решение
        SolutionState solution;
        solution = calculate_solution(ref_mg);

        logger.log("Is solved: " + std::to_string(solution.is_solved) + ' ' + solution.failure_type);

        if (solution.is_solved || solution.failure_type == 'S') {
            // logger.log("Solution norm: " + std::to_string(abs(solution.position)));
            auto ts = get_true_state_iterator(raw_mg_pas.time);
            auto error = solution.position - ts->position;
            logger.log("Error norm: " + std::to_string(abs(error)));
            if (abs(error) > 10) logger.log("HIGH ERROR");
            // auto ss_pas = pas.get_solution_state_iterator(raw_mg_pas.time);
            // auto ts_pas = pas.get_true_state_iterator(raw_mg_pas.time);
            // double error_pas = abs(ss_pas->position - ts_pas->position);
            // auto ss_act = act.get_solution_state_iterator(raw_mg_pas.time);
            // auto ts_act = act.get_true_state_iterator(raw_mg_pas.time);
            // double error_act = abs(ss_act->position - ts_act->position);
            // logger.log("Passive error: " + std::to_string(error_pas));
            // logger.log("Active error: " + std::to_string(error_act));
        }

        logger.log("");

        // low

        solution_states.push_back(solution);
    }
   
    error_type = '0';

    logger.log("Finished solving");
}

bool SatNavRel::check_raw(const RawMeasurement& raw_m) {
    if (!raw_m.is_present) {
        return false;
    }

    if (raw_m.qualflg != 0) {
        return false;
    }

    if (error_type != 'S') {
        double L1_CN0 = 20 * log10(raw_m.L1_SNR) - CN0_constant;
        double L2_CN0 = 20 * log10(raw_m.L2_SNR) - CN0_constant;
        if (L1_CN0 < CN0_min || L1_CN0 > CN0_max || L2_CN0 < CN0_min || L2_CN0 > CN0_max) {
            return false;
        }
    }

    return true;
}

RefinedMeasurement SatNavRel::refine_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act) {
    double pseudorange_delta = raw_m_pas.L1_range - raw_m_act.L1_range;
    double carrier_phase_delta = raw_m_pas.L1_phase - raw_m_act.L1_phase;

    double clock_error = pas.handler.get_clock_error(raw_m_pas.prn_id, raw_m_pas.time);
    double delay = raw_m_pas.L1_range / c + clock_error;

    GPSState gs = pas.handler.get_state(raw_m_pas.prn_id, raw_m_pas.time - delay);
    
    return {true, raw_m_pas.time, raw_m_pas.prn_id, pseudorange_delta, carrier_phase_delta, gs.position, gs.velocity}; // Минус необходим???
}

SolutionState SatNavRel::calculate_solution(const RefinedMeasurementGroupped& ref_mg) {

    SolutionState solution;
    solution.time = ref_mg.time;

    // То, что надо пронести на следующую итерацию
    DynamicFilterState dfs;
    dfs.time = ref_mg.time;

    // Находим решение пассивного КА
    auto ss_pas_it = pas.get_solution_state_iterator(ref_mg.time);
    if (ss_pas_it == pas.get_solution_states().end() || !(ss_pas_it->is_solved)) {
        solution.is_solved = false;
        solution.failure_type = 'P'; // Решение пассивного не нашлось

        dfs.is_failed = true;
        dfs_prev = dfs;

        return solution;
    }
    SolutionState ss_pas = *ss_pas_it;

    dfs.pas_pos = ss_pas.position; // Надо хранить

    // Находим решение активного КА
    auto ss_act_it = act.get_solution_state_iterator(ref_mg.time);
    if (ss_act_it == act.get_solution_states().end() || !(ss_act_it->is_solved)) {
        solution.is_solved = false;
        solution.failure_type = 'A'; // Решение активного не нашлось

        dfs.is_failed = true;
        dfs_prev = dfs;

        return solution;
    }
    SolutionState ss_act = *ss_act_it;

    // ПЕРВАЯ ИТЕРАЦИЯ И ОТСТУТСТВИЕ ОДИНОЧНЫХ РЕШЕНИЙ

    // Грубое решение для дельт
    std::vector<double> coarse = ss_act.position - ss_pas.position;

    // Проходим по всем присутствующим НКА
    for (const auto& ref_m : ref_mg.refined_measurements) {
        if (!ref_m.is_present) continue;

        std::vector<double> X_gps_minus_X_pas = ref_m.gps_position - ss_pas.position;
        double X_gps_minus_X_pas_norm = abs(X_gps_minus_X_pas);

        double delta = ref_m.gps_velocity * coarse / (c * c) + 
                       (X_gps_minus_X_pas * coarse) * (X_gps_minus_X_pas * coarse) / 
                       (X_gps_minus_X_pas_norm * X_gps_minus_X_pas_norm * X_gps_minus_X_pas_norm) / 2 -
                       coarse * coarse / X_gps_minus_X_pas_norm / 2;
        
        size_t prn_index = ref_m.prn_id - 1;
        dfs.mask[prn_index] = true;
        dfs.xi_m_1[prn_index] = ref_m.pseudorange - delta;
        dfs.xi_m_2[prn_index] = ref_m.carrier_phase;
        dfs.C_rows[prn_index] = X_gps_minus_X_pas / X_gps_minus_X_pas_norm;
    }

    // Тут завершаем первую итерацию
    if (dfs.is_first) {        
        dfs_prev = dfs;
 
        solution.is_solved = false;
        solution.failure_type = 'F';

        return solution;
    }

    // Строим компоненты вектора измерений и строки матрицы C
    size_t n = 0;
    std::vector<double> xi_m_1_temp;
    std::vector<double> d_xi_m_2_temp;
    std::vector<std::vector<double>> C_rows_temp;
    std::vector<std::vector<double>> C_rows_prev_temp;
    for (size_t prn_index = 0; prn_index < 32; prn_index++) {
        if (dfs.mask[prn_index] && dfs_prev.mask[prn_index]) {
            n++;
            xi_m_1_temp.push_back(dfs.xi_m_1[prn_index]);
            d_xi_m_2_temp.push_back(dfs.xi_m_2[prn_index] - dfs_prev.xi_m_2[prn_index]);
            C_rows_temp.push_back(dfs.C_rows[prn_index]);
            C_rows_prev_temp.push_back(dfs_prev.C_rows[prn_index]);
        }
    }

    // Собираем вектор измерений и матрицы C, C_prev
    std::vector<double> xi_m(n * 2);
    Matrix C(n, 3);
    Matrix C_prev(n, 3); // Предыдующую матрицу C надо строить заново, потому что спутники могут отличаться
    for (size_t i = 0; i < n; i++) {
        size_t j = i < n - 1 ? i + 1 : 0;
        xi_m[i] = xi_m_1_temp[i] - xi_m_1_temp[j];
        xi_m[i + n] = d_xi_m_2_temp[i] - d_xi_m_2_temp[j];
        C.at(i) = C_rows_temp[i] - C_rows_temp[j];
        C_prev.at(i) = C_rows_prev_temp[i] - C_rows_prev_temp[j];
    }

    // Строим матрицу C0
    Matrix C0(n * 2, 6, 0.0);
    C0.insert(C, 0, 0);
    C0.insert(C, n, 3);

    // Нужные вещи
    double dt = dfs.time - dfs_prev.time;
    Matrix E = identity(3);
    Matrix Omega = {{{0, earth_rotation_rate, 0}, {-earth_rotation_rate, 0, 0}, {0, 0, 0}}};
    double X_norm = abs(dfs_prev.pas_pos);
    double omega2 = earth_mu / (X_norm * X_norm * X_norm);
    Matrix X_prod_X_T = col_by_row(dfs_prev.pas_pos, dfs_prev.pas_pos);
    logger.log(std::to_string(abs(dfs.pas_pos - dfs_prev.pas_pos)));
    Matrix derivative = (Omega * Omega + (E - X_prod_X_T / (X_norm * X_norm)  * 3) * omega2) * dt * dt * (-1);

    // Строим матрицу, обратную матрице B
    Matrix B1(6, 6);
    B1.insert(E, 0, 0);
    B1.insert(E * (-1), 0, 3);
    B1.insert(derivative * (-1), 3, 0);
    B1.insert(E + derivative - Omega * dt * 2, 3, 3);

    logger.log("dt: " + std::to_string(dt));
    if (dt > 10) std::cout << dfs.time << " -- " << dt << std::endl;

    // Строим оценку вектора состояния путем прогноза
    double d = dot(dfs_prev.pas_pos, dfs_prev.x);
    double r2 = dot(dfs_prev.pas_pos, dfs_prev.pas_pos);
    double del = 3 * d / r2; // + 3 / 2 * d * d / (r2 * r2) + 3 / 2 * dot(dfs_prev.x, dfs_prev.x) / (r2 * r2);

    std::vector<double> dx_est = (E + Omega * dt * 2) * dfs_prev.dx - Omega * Omega * dfs_prev.x * dt * dt - (dfs_prev.x - dfs_prev.pas_pos * del) * omega2 * dt * dt;
    std::vector<double> x_est = dfs_prev.x + dx_est;

    // Строим оценку вектора измерений (3.55)
    std::vector<double> x_m_est = C * x_est;
    std::vector<double> dx_m_est = (C - C_prev) * dfs_prev.x + C * dx_est;

    // Собираем все в векторы кси     
    std::vector<double> xi_est(6);
    for (size_t i = 0; i < 3; i++) {
        xi_est[i] = x_est[i];
        xi_est[i + 3] = dx_est[i];
    }
    std::vector<double> xi_m_est(n * 2);
    for (size_t i = 0; i < n; i++) {
        xi_m_est[i] = x_m_est[i];
        xi_m_est[i + n] = dx_m_est[i];
    }

    // Сглаживание измерений
    for (size_t i = 0; i < n; i++) {
        xi_m[i] = 1 / T_p * xi_m[i] + (1 - 1 / T_p) * xi_m[i + n];
        xi_m_est[i] = 1 / T_p * xi_m_est[i] + (1 - 1 / T_p) * xi_m_est[i + n];
    }

    // Динамическая фильтрация
    std::vector<double> P = C0.transpose() * (xi_m - xi_m_est);
    W = B1.transpose() * lambda * W * lambda * B1 + C0.transpose() * C0;

    // Проверка обратимости W
    std::vector<double> d_xi(6);
    try {
        d_xi = W.inverse() * P;
    } catch (std::runtime_error) {
        logger.log("W необратима");
    }
    std::vector<double> xi = xi_est + d_xi;

    // Сохраняем на следующую итерацию
    for (size_t i = 0; i < 3; i++) {
        dfs.x[i] = xi[i];
        dfs.dx[i] = xi[i + 3];
    }
    dfs_prev = dfs;

    solution.position = dfs.x;
    solution.is_solved = true;

    return solution;
}

const std::vector<State>& SatNavRel::get_true_states() const {
    return true_states;
}

const std::vector<SolutionState>& SatNavRel::get_solution_states() const {
    return solution_states;
}

std::vector<State>::const_iterator SatNavRel::get_true_state_iterator(int time) const {
    size_t n = true_states.size();

    size_t lo = 0, hi = n - 1;
    while (lo <= hi) {
        size_t mid = lo + (hi - lo) / 2;

        if (true_states[mid].time == time) {
            return true_states.cbegin() + mid;
        }

        if (true_states[mid].time < time) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    return true_states.cend();
}
