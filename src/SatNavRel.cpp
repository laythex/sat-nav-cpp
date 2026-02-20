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
            if (!pas.check_raw(raw_mg_pas.raw_measurements[prn_index]) || !act.check_raw(raw_mg_act.raw_measurements[prn_index])) continue;

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
            // logger.log("Solution x: " + std::to_string(solution.position[0]));
            // logger.log("Solution y: " + std::to_string(solution.position[1]));
            // logger.log("Solution z: " + std::to_string(solution.position[2]));

            auto ts = get_true_state_iterator(raw_mg_pas.time);
            // logger.log("True norm: " + std::to_string(abs(ts->position)));
            // logger.log("True x: " + std::to_string(ts->position[0]));
            // logger.log("True y: " + std::to_string(ts->position[1]));
            // logger.log("True z: " + std::to_string(ts->position[2]));

            auto error = solution.position - ts->position;
            logger.log("Error norm: " + std::to_string(abs(error)));
            // logger.log("Error x: " + std::to_string(error[0]));
            // logger.log("Error y: " + std::to_string(error[1]));
            // logger.log("Error z: " + std::to_string(error[2]));
            if (abs(error) > 15) logger.log("HIGH ERROR");
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
    if (dfs_prev.time == 0) dfs.is_first = true; // Самая первая итерация
    if (dfs_prev.is_failed) dfs.is_first = true; // Последняя итерация не удалась - начинаем заново
    if (dfs_prev.is_first) dfs.is_second = true; // полукостыль

    // Находим решение пассивного КА
    auto ss_pas_it = pas.get_solution_state_iterator(ref_mg.time);
    if (ss_pas_it == pas.get_solution_states().end() || !(ss_pas_it->is_solved)) { // Лишние скобки?
        solution.is_solved = false;
        solution.failure_type = 'P'; // Решение пассивного не нашлось

        dfs.is_failed = true;
        dfs_prev = dfs;

        return solution;
    }
    SolutionState ss_pas = *ss_pas_it;

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

    // Грубое решение для дельт
    std::vector<double> coarse = ss_act.position - ss_pas.position;

    // То, что не надо проносить на следующую итерацию
    std::vector<double> xi_m_1_full(32);
    std::vector<std::vector<double>> C_rows_full(32);

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
        dfs.xi_m_2[prn_index] = ref_m.carrier_phase;
        xi_m_1_full[prn_index] = ref_m.pseudorange - delta;
        C_rows_full[prn_index] = X_gps_minus_X_pas / X_gps_minus_X_pas_norm;
    }

    /*
        Предпредыдущий шаг:
        xi_m_2

        Предыдущий шаг:
        xi_m_1, d_xi_m_2, C => x, dx - из системы

        Текущий шаг:
        x_m_1, d_xi_m_2, С
        из прогноза x, dx
        решаем ДФ
    */

    // Обрабатываем предпредыдущий шаг
    if (dfs.is_first) {
        dfs_prev = dfs;

        // К этому моменту в дфс лежат маска и xi_m_2
        solution.is_solved = false;
        solution.position = dfs.x; // для отладки (для отображения в логах)
        solution.failure_type = 'F';
        return solution;
    }

    // Строим компоненты вектора измерений и строки матрицы C
    size_t n = 0;
    std::vector<double> xi_m_1;
    std::vector<double> d_xi_m_2;
    std::vector<std::vector<double>> C_rows;
    for (size_t prn_index = 0; prn_index < 32; prn_index++) {
        if (dfs.mask[prn_index] && dfs_prev.mask[prn_index]) {
            n++;
            xi_m_1.push_back(xi_m_1_full[prn_index]);
            d_xi_m_2.push_back(dfs.xi_m_2[prn_index] - dfs_prev.xi_m_2[prn_index]);
            C_rows.push_back(C_rows_full[prn_index]);
        }
    }

    // Собираем вектор измерений и матрицу C
    std::vector<double> xi_m(n * 2);
    Matrix C(n, 3);
    for (size_t i = 0; i < n; i++) {
        size_t j = i < n - 1 ? i + 1 : 0;
        xi_m[i] = xi_m_1[i] - xi_m_1[j];
        xi_m[i + n] = d_xi_m_2[i] - d_xi_m_2[j];
        C.at(i) = C_rows[i] - C_rows[j];
    }

    // Матрица C может оказаться плохой 
    try {
        Matrix C1 = (C.transpose() * C).inverse();
        double c = sqrt(C1.trace());
        if (c > 5) {
            dfs.is_failed = true;
            dfs_prev = dfs;

            solution.is_solved = false;
            solution.failure_type = 'c';
            return solution;  
        }
    } catch (std::runtime_error e) {
        dfs.is_failed = true;
        dfs_prev = dfs;

        solution.is_solved = false;
        solution.failure_type = 'C';
        return solution; 
    }

    // Строим матрицу C0
    Matrix C0(n * 2, 6, 0.0);
    C0.insert(C, 0, 0);
    C0.insert(C, n, 3);

    // Обрабатываем предыдущий шаг
    if (dfs.is_second) {
        std::vector<double> xi = (C.transpose() * C).inverse() * C.transpose() * xi_m; // На этом шаге пока без прогнозов и ДФ
        for (size_t i = 0; i < 3; i++) {
            dfs.x[i] = xi[i];
            dfs.dx[i] = xi[i + 3];
        }
        dfs_prev = dfs;
        W = Matrix(6, 6, 0.0); // Обнуляем матрицу W

        // К этому моменту в дфс лежит маска, xi_m_2, x и dx => готовы начать со следующего шага
        solution.is_solved = false;
        solution.position = dfs.x; // для отладки (для отображения в логах)
        solution.failure_type = 'S';
        return solution;
    }

    // Нужные вещи
    double dt = dfs.time - dfs_prev.time;
    Matrix E = identity(3);
    Matrix Omega = {{{0, earth_rotation_rate, 0}, {-earth_rotation_rate, 0, 0}, {0, 0, 0}}}; // Сделать нормальную матрицу поворота?
    double X_norm = abs(ss_pas.position);
    double omega2 = earth_mu / (X_norm * X_norm * X_norm);
    Matrix X_prod_X_T = col_by_row(ss_pas.position, ss_pas.position);
    Matrix derivative = (Omega * Omega + (E - X_prod_X_T / (X_norm * X_norm)  * 3) * omega2) * dt * dt * (-1);

    logger.log("dt: " + std::to_string(dt));
    if (dt > 10) std::cout << dfs.time << std::endl;

    // Строим оценку вектора состояния путем прогноза
    std::vector<double> dx_est = (E + Omega * dt * 2) * dfs_prev.dx - Omega * Omega * dfs_prev.x * dt * dt; 
                                //  (E - X_prod_X_T / (X_norm * X_norm) * 3) * dfs_prev.x * omega2 * dt * dt;
    std::vector<double> x_est = dfs_prev.x + dx_est;

    // Чистое интегрирование в приращениях
    dfs.x = x_est;
    dfs.dx = dx_est;
    dfs_prev = dfs;
    solution.position = x_est;
    solution.is_solved = true;
    return solution;

    // Строим оценку вектора состояния
    std::vector<double> xi_est(6);
    for (size_t i = 0; i < 3; i++) {
        xi_est[i] = x_est[i];
        xi_est[i + 3] = dx_est[i];
    }

    // Строим оценку вектора измерений из оценки вектора состояния
    std::vector<double> xi_m_est = C0 * xi_est;

    logger.log("prev dx: " + std::to_string(abs(dfs_prev.dx)));
    logger.log("prev x:  " + std::to_string(abs(dfs_prev.x)));             
    logger.log("dx_est:  " + std::to_string(abs(dx_est)));
    logger.log("d_xi_m:  " + std::to_string(abs(xi_m - xi_m_est))); 

    // Строим матрицу, обратную матрице B
    Matrix B1(6, 6);
    B1.insert(E, 0, 0);
    B1.insert(E * (-1), 0, 3);
    B1.insert(derivative * (-1), 3, 0);
    B1.insert(E + derivative - Omega * dt * 2, 3, 3);

    // Динамическая фильтрация
    std::vector<double> P = C0.transpose() * (xi_m - xi_m_est);
    W = B1.transpose() * lambda * W * lambda * B1 + C0.transpose() * C0;
    std::vector<double> d_xi = W.inverse() * P;
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
