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

        if (solution.is_solved) {
            logger.log("GDOP: " + std::to_string(solution.GDOP));
            logger.log("Solution norm: " + std::to_string(abs(solution.position)));

            auto ts = get_true_state_iterator(raw_mg_pas.time);
            double error = abs(solution.position - ts->position);
            logger.log("Error: " + std::to_string(error));
            if (error > 5) logger.log("HIGH ERROR");
            auto ss_pas = pas.get_solution_state_iterator(raw_mg_pas.time);
            auto ts_pas = pas.get_true_state_iterator(raw_mg_pas.time);
            double error_pas = abs(ss_pas->position - ts_pas->position);
            auto ss_act = act.get_solution_state_iterator(raw_mg_pas.time);
            auto ts_act = act.get_true_state_iterator(raw_mg_pas.time);
            double error_act = abs(ss_act->position - ts_act->position);
            logger.log("Passive error: " + std::to_string(error_pas));
            logger.log("Active error: " + std::to_string(error_act));
        }

        logger.log("");

        // low

        solution_states.push_back(solution);
    }
   
    error_type = '0';

    logger.log("Finished solving");
}

RefinedMeasurement SatNavRel::refine_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act) {
    double pseudorange_delta = raw_m_act.L1_range - raw_m_pas.L1_range;
    double carrier_phase_delta = raw_m_act.L1_phase - raw_m_pas.L1_phase;

    double clock_error = pas.handler.get_clock_error(raw_m_pas.prn_id, raw_m_pas.time);
    double delay = raw_m_pas.L1_range / c + clock_error;

    GPSState gs = pas.handler.get_state(raw_m_pas.prn_id, raw_m_pas.time - delay);
    
    return {true, raw_m_pas.time, raw_m_pas.prn_id, pseudorange_delta, carrier_phase_delta, gs.position, gs.velocity};
}

SolutionState SatNavRel::calculate_solution(const RefinedMeasurementGroupped& ref_mg) {

    SolutionState solution;
    solution.time = ref_mg.time;

    // Находим решение пассивного КА
    auto ss_pas_it = pas.get_solution_state_iterator(ref_mg.time);
    if (ss_pas_it == pas.get_solution_states().end() || !(ss_pas_it->is_solved)) { // Лишние скобки?
        solution.is_solved = false;
        solution.failure_type = 'P';
        return solution;    // Решение не нашлось
    }
    SolutionState ss_pas = *ss_pas_it;

    // Находим решение активного КА
    auto ss_act_it = act.get_solution_state_iterator(ref_mg.time);
    if (ss_act_it == act.get_solution_states().end() || !(ss_act_it->is_solved)) {
        solution.is_solved = false;
        solution.failure_type = 'A';
        return solution;    // Решение не нашлось
    }
    SolutionState ss_act = *ss_act_it;

    // Грубое решение
    std::vector<double> coarse = ss_act.position - ss_pas.position;

    // Проходим по всем присутствующим НКА
    size_t n = 0;
    std::vector<double> pr_act_minus_pr_pas_minus_delta;
    std::vector<double> cp_act_minus_cp_pas;
    std::vector<std::vector<double>> gps_pos_minus_pas_pos_norm;
    std::vector<double> pr_act_minus_pr_pas_minus_delta_est;
    std::vector<double> cp_act_minus_cp_pas_est;
    for (const auto& ref_m : ref_mg.refined_measurements) {
        if (!ref_m.is_present) continue;

        std::vector<double> X_gps_minus_X_act = ref_m.gps_position - ss_act.position;
        std::vector<double> X_gps_minus_X_pas = ref_m.gps_position - ss_pas.position;

        double delta = ref_m.gps_velocity * coarse / (c * c) + 
                       (X_gps_minus_X_pas * coarse) * (X_gps_minus_X_pas * coarse) / 
                       (abs(X_gps_minus_X_pas) * abs(X_gps_minus_X_pas) * abs(X_gps_minus_X_pas)) / 2 -
                       coarse * coarse / abs(X_gps_minus_X_pas) / 2;

        n++;
        
        // Из точных решений
        pr_act_minus_pr_pas_minus_delta.push_back(ref_m.pseudorange - delta);
        cp_act_minus_cp_pas.push_back(ref_m.carrier_phase);
        gps_pos_minus_pas_pos_norm.push_back(X_gps_minus_X_pas / abs(X_gps_minus_X_pas));

        // Из приближений КСВ
        pr_act_minus_pr_pas_minus_delta_est.push_back(abs(X_gps_minus_X_act) - abs(X_gps_minus_X_pas) - delta);
        cp_act_minus_cp_pas_est.push_back(abs(X_gps_minus_X_act) - abs(X_gps_minus_X_pas));
    }

    // Строим векторы измерений - истинный и приблизительный
    std::vector<double> xi_m_1(n);
    std::vector<double> xi_m_2(n);
    std::vector<double> xi_m_1_est(n); // почему n?
    std::vector<double> xi_m_2_est(n);
    for (size_t i = 0; i < n; i++) {
        size_t j = i < n - 1 ? i + 1 : 0;

        xi_m_1[i] = pr_act_minus_pr_pas_minus_delta[i] - pr_act_minus_pr_pas_minus_delta[j];
        xi_m_2[i] = cp_act_minus_cp_pas[i] - cp_act_minus_cp_pas[j];

        xi_m_1_est[i] = pr_act_minus_pr_pas_minus_delta_est[i] - pr_act_minus_pr_pas_minus_delta_est[j];
        xi_m_2_est[i] = cp_act_minus_cp_pas_est[i] - cp_act_minus_cp_pas_est[j];
    }

    // Обрабатываем первое решение
    if (time_prev == 0) {
        solution.position = coarse;
        solution.is_solved = true;

        time_prev = ref_mg.time;
        coarse_prev = coarse;
        xi_m_2_prev = xi_m_2;
        xi_m_2_est_prev = xi_m_2_est;

        return solution;
    }

    // Достраиваем вектор измерений
    std::vector<double> xi_m(n * 2);
    std::vector<double> xi_m_est(n * 2);
    for (size_t i = 0; i < n; i++) {
        xi_m[i] = xi_m_1[i];
        xi_m[i + n] = xi_m_2[i] - xi_m_2_prev[i];

        xi_m_est[i] = xi_m_1_est[i];
        xi_m_est[i + n] = xi_m_2_est[i] - xi_m_2_est_prev[i];
    }

    // Строим приблизительный вектор состояния 
    std::vector<double> xi_est(6);
    for (size_t i = 0; i < 3; i++) {
        xi_est[i] = coarse[i];
        xi_est[i + 3] = coarse[i] - coarse_prev[i];
    }

    // Строим матрицу C0
    Matrix C0(n * 2, 6, 0.0);
    Matrix C(n, 3);
    for (size_t i = 0; i < n; i++) {
        size_t j = i < n - 1 ? i + 1 : 0;
        C.at(i) = gps_pos_minus_pas_pos_norm[i] - gps_pos_minus_pas_pos_norm[j];
    }
    C0.insert(C, 0, 0);
    C0.insert(C, 3, 3);

    double dt = ref_mg.time - time_prev;

    Matrix E = identity(3);
    double omega = earth_rotation_rate;
    Matrix Omega = {{{0, omega, 0}, {omega, 0, 0}, {0, 0, 0}}};

    double X_norm2 = abs(ss_pas.position) * abs(ss_pas.position);
    Matrix X_prod_X_T = col_by_row(ss_pas.position, ss_pas.position);
    Matrix derivative = (Omega * Omega + (E - X_prod_X_T / X_norm2 * 3) * omega * omega) * dt * dt * (-1);

    // Строим матрицу B
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

    for (size_t i = 0; i < 3; i++) {
        solution.position[i] = xi[i];
    }
    solution.is_solved = true;

    // Сохраняем на следующую итерацию
    time_prev = ref_mg.time;
    coarse_prev = coarse;
    xi_m_2_prev = xi_m_2;
    xi_m_2_est_prev = xi_m_2_est;

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
