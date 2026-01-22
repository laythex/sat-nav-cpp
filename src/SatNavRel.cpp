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

    // Заполняем вектор грубых решений // УБРАТЬ?
    for (const SolutionState& ss_pas : pas.get_solution_states()) {

        auto ss_act_it = act.get_solution_state_iterator(ss_pas.time);
        if (ss_act_it == act.get_solution_states().end()) continue;
        SolutionState ss_act = *ss_act_it;
        
        if (ss_pas.is_solved && ss_act.is_solved) {
            coarse_solution_states.push_back({ss_pas.time, ss_act.position - ss_pas.position, ss_act.velocity - ss_pas.velocity, true, '0', 0});
        } else {
            coarse_solution_states.push_back({ss_pas.time, {0, 0, 0}, {0, 0, 0}, false, '0', 0});
        }
    }
}

SatNavRel::SatNavRel(const SatNavRel& sn) : pas(sn.pas), act(sn.act),
                                            true_states(sn.true_states), 
                                            coarse_solution_states(sn.coarse_solution_states) {}

void SatNavRel::solve_relative(char et, int ti, int tf) {

    logger.log("Beginning to solve...");

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

        // Находим грубую оценку относительного вектора из КСВ
        auto coarse_solution_it = get_coarse_solution_state_iterator(raw_mg_pas.time); 

        // Начинаем динамическую фильтрацию
        double lambda_x = T_x / (T_x + 1);
        double lambda_v = T_v / (T_v + 1);
        lambda.insert(identity(3) * lambda_x, 0, 0);
        lambda.insert(identity(3) * lambda_v, 3, 3);

        // Ищем решение
        SolutionState solution;
        // if (coarse_solution_it != coarse_solution_states.end() && coarse_solution_it->is_solved) { // Если грубое решение нашлось
        //     solution = calculate_solution(ref_mg, coarse_solution_it->position);
        // } else { // Иначе решаем два раза для уточнения
        //     solution = calculate_solution(ref_mg, {0, 0, 0});
        //     solution = calculate_solution(ref_mg, solution.position);
        // }

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

        // if (error_type != 'L') {
        //     if (solution.is_solved) {
        //         std::vector<unsigned> low = check_low(solution, ref_mg);
        //         if (low.size() > 0) {
        //             for (unsigned prn_id : low) {
        //                 unsigned prn_index = prn_id - 1;
        //                 ref_mg.refined_measurements[prn_index].is_present = false;
        //             }
        //             solution = calculate_solution(ref_mg);
        //         }
        //     }
        // }

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

    // Находим последнее решение
    

    // Проходим по всем присутствующим НКА
    size_t n = 0;
    std::vector<double> pr_act_minus_pr_pas_minus_delta;
    std::vector<double> cp_act_minus_cp_pas;
    std::vector<std::vector<double>> gps_pos_minus_pas_pos;
    std::vector<std::vector<double>> gps_pos_minus_pas_pos_norm;
    for (const auto& ref_m : ref_mg.refined_measurements) {
        if (ref_m.is_present) {
            std::vector<double> X_gps_m_X_pas = ref_m.gps_position - ss_pas.position;

            double delta = ref_m.gps_velocity * coarse / (c * c) + 
                           (X_gps_m_X_pas * coarse) * (X_gps_m_X_pas * coarse) / (abs(X_gps_m_X_pas) * abs(X_gps_m_X_pas) * abs(X_gps_m_X_pas)) / 2 -
                           coarse * coarse / abs(X_gps_m_X_pas) / 2;

            n++;
            pr_act_minus_pr_pas_minus_delta.push_back(ref_m.pseudorange - delta);
            cp_act_minus_cp_pas.push_back(ref_m.carrier_phase);
            gps_pos_minus_pas_pos.push_back(X_gps_m_X_pas);
            gps_pos_minus_pas_pos_norm.push_back(X_gps_m_X_pas / abs(X_gps_m_X_pas));
        }
    }

    // Строим необходимые векторы и матрицы
    std::vector<double> xi_m(n * 2);
    Matrix C0(n * 2, 6, 0.0);
    Matrix B(6, 6);

    double dt; // НАДО

    Matrix C(n, 3);
    for (size_t i = 0; i < n; i++) {
        size_t j = i < n - 1 ? i + 1 : 0;

        xi_m[i] = pr_act_minus_pr_pas_minus_delta[i] - pr_act_minus_pr_pas_minus_delta[j];
        xi_m[i + n] = cp_act_minus_cp_pas[i] - cp_act_minus_cp_pas[j];

        C.at(i) = gps_pos_minus_pas_pos_norm[i] - gps_pos_minus_pas_pos_norm[j];
    }

    C0.insert(C, 0, 0);
    C0.insert(C, 3, 3);

    Matrix E = identity(3);
    double omega = earth_rotation_rate;
    Matrix Omega = {{{0, omega, 0}, {omega, 0, 0}, {0, 0, 0}}};

    double X_norm2 = abs(ss_pas.position) * abs(ss_pas.position);
    Matrix X_prod_X_T = col_by_row(ss_pas.position, ss_pas.position);
    Matrix derivative = (Omega * Omega + (E - X_prod_X_T / X_norm2 * 3) * omega * omega) * dt * dt * (-1);

    B.insert(E, 0, 0);
    B.insert(E * (-1), 0, 3);
    B.insert(derivative * (-1), 3, 0);
    B.insert(E + derivative - Omega * dt * 2, 3, 3);

    std::vector<double> xi_est = estimate_state(ss_pas, ss_last);
    std::vector<double> xi_m_est = estimate_measurement();

    std::vector<double> P = C0.transpose() * (xi_m - xi_m_est);
    W = B.transpose() * lambda * W * lambda * B + C0.transpose() * C0;
    std::vector<double> d_xi = W.inverse() * P;
    std::vector<double> xi = xi_est + d_xi;

    // solution.position = dX;
    solution.is_solved = true;

    return solution;
}

State SatNavRel::estimate_state(const SolutionState& ss_pas, const SolutionState& ss_last, int time) {

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

std::vector<SolutionState>::const_iterator SatNavRel::get_coarse_solution_state_iterator(int time) const {
    size_t n = coarse_solution_states.size();

    size_t lo = 0, hi = n - 1;
    while (lo <= hi) {
        size_t mid = lo + (hi - lo) / 2;

        if (coarse_solution_states[mid].time == time) {
            return coarse_solution_states.cbegin() + mid;
        }

        if (coarse_solution_states[mid].time < time) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    return coarse_solution_states.cend();
}
