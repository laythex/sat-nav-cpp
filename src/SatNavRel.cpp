#include "SatNavRel.hpp"
#include "SatNav.hpp"

SatNavRel::SatNavRel(SatNav& passive, SatNav& active) : pas(passive), act(active), logger("RGPS") {}

SatNavRel::SatNavRel(const SatNavRel& sn) : pas(sn.pas), act(sn.act), true_states(sn.true_states), logger("RGPS") {}
                                        
void SatNavRel::solve_separately(char et, unsigned ti, unsigned tf) {
    pas.solve(et, ti, tf);
    act.solve(et, ti, tf);

    for (const State& ts_pas : pas.get_true_states()) {
        unsigned time = ts_pas.time;
        auto ts_act_it = act.get_true_state_iterator(time);
        if (ts_act_it == act.get_true_states().end()) continue;
        State ts_act = *ts_act_it;
        true_states.push_back({time, ts_act.position - ts_pas.position, ts_act.velocity - ts_pas.velocity});
    }
}

void SatNavRel::solve_relative(char et, unsigned ti, unsigned tf) {

    solve_separately(et, ti, tf);

    logger.log("Starting RGPS solving");

    error_type = et;
    double err_avg = 0;
    double err_cnt = 0;

    unsigned t_front = std::min(pas.raw_measurements_groupped.front().time,
                                act.raw_measurements_groupped.front().time);
    unsigned t_back = std::max(pas.raw_measurements_groupped.back().time,
                               act.raw_measurements_groupped.back().time);

    for (unsigned time = t_front; time < t_back; time += dt) {

        if (ti > 0 || tf > 0) {
            if ((time - t_front) < ti || (time - t_front) >= tf) continue;
        }

        RawMeasurementGroupped raw_mg_pas = {time, std::vector<RawMeasurement>(32)};
        auto raw_mg_pas_it = pas.get_raw_measurement_groupped_iterator(time);
        if (raw_mg_pas_it != pas.raw_measurements_groupped.end()) {
            raw_mg_pas = *raw_mg_pas_it;
        }

        RawMeasurementGroupped raw_mg_act = {time, std::vector<RawMeasurement>(32)};
        auto raw_mg_act_it = act.get_raw_measurement_groupped_iterator(time);
        if (raw_mg_act_it != act.raw_measurements_groupped.end()) {
            raw_mg_act = *raw_mg_act_it;
        }

        logger.logv("Time", time);

        // Обрабатываем сырые измерения с каждого доступного спутника
        std::vector<RefinedMeasurement> ref_ms(32);
        for (unsigned prn_id = 1; prn_id <= 32; ++prn_id) {
            unsigned prn_index = prn_id - 1;

            RawMeasurement raw_m_pas = raw_mg_pas.raw_measurements[prn_index];
            RawMeasurement raw_m_act = raw_mg_act.raw_measurements[prn_index];

            // Проверяем, что с измерениями все ок
            if (!check_raw(raw_m_pas, raw_m_act)) continue;

            logger.logv("Refining", prn_id);

            // Непосредственно их обрабатываем
            ref_ms[prn_index] = refine_raw(raw_m_pas, raw_m_act);
        }

        // Группируем обработанные измерения со всех спутников
        RefinedMeasurementGroupped ref_mg = {time, ref_ms};
        refined_measurements_groupped.push_back(ref_mg);

        // Ищем решение
        SolutionState solution = calculate_solution(ref_mg);
        solution_states.push_back(solution);

        logger.log("Is solved: " + std::to_string(solution.is_solved) + ' ' + solution.failure_type);
        auto ts = get_true_state_iterator(raw_mg_pas.time);
        double error = abs(solution.position - ts->position);

        logger.logv("Error norm", error);
        if (solution.is_solved) {
            err_avg += error;
            ++err_cnt;
            if (error > 3) {
                logger.log("HIGH ERROR");
            }
        }
        logger.lnbr();
    }
   
    error_type = '0';

    err_avg /= static_cast<double>(err_cnt);
    logger.logv("Average error", err_avg);
    std::cout << std::format("Avergate error: {:.3f}, Solution count: {}", err_avg, err_cnt) << std::endl;

    logger.log("Finished RGPS solving");
}

bool SatNavRel::check_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act) {
    if (!raw_m_pas.is_present || !raw_m_act.is_present) {
        return false;
    }

    if (raw_m_pas.qualflg != 0 || raw_m_act.qualflg != 0) {
        return false;
    }

    if (error_type != 'S') {
        double min_CN0 = std::min(raw_m_pas.L1_CN0, raw_m_act.L1_CN0);
        double max_CN0 = std::max(raw_m_pas.L1_CN0, raw_m_act.L1_CN0);

        if (min_CN0 < CN0_min_threshold || max_CN0 > CN0_max_threshold) {
            logger.log("Bad CN0: " + std::to_string(raw_m_pas.prn_id) + ", " + std::to_string(min_CN0) + ", " + std::to_string(max_CN0));
            return false;
        }
    }

    if (error_type != 'F') {
        if (pas.check_fading(raw_m_pas) || act.check_fading(raw_m_act)) {
            logger.logv("Fading", raw_m_pas.prn_id);
            return false;
        }
    }

    return true;
}

RefinedMeasurement SatNavRel::refine_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act) {
    double pseudorange_delta = raw_m_pas.L1_range - raw_m_act.L1_range;
    double carrier_phase_delta = raw_m_pas.L1_phase - raw_m_act.L1_phase;

    double delay = raw_m_pas.L1_range / c;
    Matrix R = (E3 + Omega * delay);

    State gs = pas.handler.get_state(raw_m_pas.prn_id, raw_m_pas.time);

    return {true, 
            raw_m_pas.time, 
            raw_m_pas.prn_id, 
            pseudorange_delta, 
            carrier_phase_delta, 
            R * gs.position, 
            R * gs.velocity};
}

SolutionState SatNavRel::calculate_solution(const RefinedMeasurementGroupped& ref_mg) {
    unsigned time = ref_mg.time;

    SolutionState solution;
    solution.time = time;

    DynamicFilterState dfs;
    dfs.time = time;
    
    std::vector<double> pas_pos = get_coarse_position(SatType::PASSIVE, time);
    std::vector<double> act_pos = get_coarse_position(SatType::ACTIVE, time);
    std::vector<double> coarse = act_pos - pas_pos;
    dfs.pas_pos = pas_pos;

    // Проходим по всем присутствующим НКА
    for (const auto& ref_m : ref_mg.refined_measurements) {

        if (!ref_m.is_present) continue;
        
        double delta = calculate_delta(ref_m.gps_position - dfs.pas_pos, coarse, ref_m.gps_velocity);

        size_t prn_index = ref_m.prn_id - 1;
        dfs.mask[prn_index] = true;
        dfs.xi_m_1[prn_index] = ref_m.pseudorange + delta;
        dfs.xi_m_2[prn_index] = ref_m.carrier_phase + delta;
        dfs.C_rows[prn_index] = normalize(ref_m.gps_position - dfs.pas_pos);
    }

    if (df_state == 0) {
        dfs.x = coarse;
        dfs_prev = dfs;
        df_state = 1;

        solution.position = dfs.x;
        solution.is_solved = false;
        solution.failure_type = '0';
        return solution;
    }

    if (df_state == 1) {
        dfs.x = coarse;
        dfs.dx = coarse - dfs_prev.x;
        dfs.x_est = coarse;
        dfs_prev = dfs;
        df_state = 2;

        solution.position = dfs.x;
        solution.is_solved = false;
        solution.failure_type = '1';
        return solution;
    }

    // Строим компоненты вектора измерений и строки матрицы C
    size_t n = 0;
    std::vector<double> xi_m_1_temp;
    std::vector<double> d_xi_m_2_temp;
    std::vector<std::vector<double>> C_rows_temp;
    std::vector<std::vector<double>> dC_rows_temp;
    logger.log("Working with:", ' ');
    for (size_t prn_id = 1; prn_id < 33; ++prn_id) {
        size_t prn_index = prn_id - 1;
        if (dfs.mask[prn_index] && dfs_prev.mask[prn_index]) {
            ++n;
            xi_m_1_temp.push_back(dfs.xi_m_1[prn_index]);
            d_xi_m_2_temp.push_back(dfs.xi_m_2[prn_index] - dfs_prev.xi_m_2[prn_index]);
            C_rows_temp.push_back(dfs.C_rows[prn_index]);
            dC_rows_temp.push_back(dfs.C_rows[prn_index] - dfs_prev.C_rows[prn_index]);
            logger.log(prn_id, ' ');
        }
    }
    logger.lnbr();
    logger.logv("Sat count (n)", n);
    
    // Собираем вектор измерений и матрицы C, dC
    std::vector<double> xi_m(n * 2);
    Matrix C(n, 3);
    Matrix dC(n, 3);
    for (size_t i = 0; i < n; ++i) {
        size_t j = i < n - 1 ? i + 1 : 0;
        xi_m[i] = xi_m_1_temp[i] - xi_m_1_temp[j];
        xi_m[i + n] = d_xi_m_2_temp[i] - d_xi_m_2_temp[j];
        C.at(i) = C_rows_temp[i] - C_rows_temp[j];
        dC.at(i) = dC_rows_temp[i] - dC_rows_temp[j];
    }

    Matrix B1 = calculate_B1();
    Matrix C0 = zero(n * 2, 6);
    C0.insert(C, 0, 0);
    C0.insert(C, n, 3);

    // Строим оценку вектора состояния путем прогноза
    std::vector<double> dx_est = estimate_dx();
    std::vector<double> x_est = dfs_prev.x + dx_est;
    dfs.x_est = x_est; // Запоминаем для следующей итерации

    logger.logv("Model error (x_est)", abs(x_est - get_true_state_iterator(time)->position));
    logger.logv("Model error (dx_est)", abs(dx_est - (get_true_state_iterator(time)->position - get_true_state_iterator(time - dt)->position)));

    // Строим оценку вектора измерений (3.55)
    std::vector<double> x_m_est = C * x_est;
    std::vector<double> dx_m_est = dC * dfs_prev.x_est + C * dx_est;

    // Собираем все в векторы кси     
    std::vector<double> xi_est(6);
    for (size_t i = 0; i < 3; ++i) {
        xi_est[i] = x_est[i];
        xi_est[i + 3] = dx_est[i];
    }

    std::vector<double> xi_m_diff(n * 2);
    for (size_t i = 0; i < n; ++i) {
        xi_m_diff[i] = xi_m[i] - x_m_est[i];
        xi_m_diff[i + n] = xi_m[i + n] - dx_m_est[i];
    }

    // std::cout << time << std::endl << xi_m << std::endl << xi_m_est << std::endl;
    // std::cout << (C0.transpose() * C0).inverse() * C0.transpose() * xi_m << std::endl;
    // std::cout << (C0.transpose() * C0).inverse() * C0.transpose() * xi_m_est << std::endl;
    // std::cout << get_true_state_iterator(time)->position << get_true_state_iterator(time)->position - get_true_state_iterator(time - dt)->position << std::endl;
    // std::cout << std::endl;

    if (is_pure_modeling_mode) {
        dfs.x = x_est;
        dfs.dx = dx_est;
        dfs_prev = dfs;
        is_model_step = true;
        ++consecutive_model_steps;
        since_last_model_step = 0;
        solution.position = dfs.x;
        solution.is_solved = true;
        return solution;
    }

    // 3.68
    while (n > 0) {
        if (consecutive_model_steps > model_steps_meas_diff_threshold) {
            logger.logv("Model steps exceed threshold, no meas diff filtering", consecutive_model_steps);
            break;
        }

        if (since_last_model_step < model_steps_relaxation_threshold) {
            logger.logv("Too early since last model step, no meas diff filtering", since_last_model_step);
            break;
        }

        size_t max_at = find_max_abs(xi_m_diff);
        double max_diff = std::abs(xi_m_diff[max_at]);

        if (max_diff > meas_diff_threshold) {
            if (max_at > n) max_at -= n;

            std::vector<double> xi_m_diff_new((n - 1) * 2);
            Matrix C_new(n - 1, 3);
            
            for (size_t i = 0; i < n - 1; ++i) {
                if (i < max_at) {
                    xi_m_diff_new[i] = xi_m_diff[i];
                    xi_m_diff_new[i + n] = xi_m_diff[i + n];
                    C_new.at(i) = C(i);
                } else {
                    xi_m_diff_new[i] = xi_m_diff[i + 1];
                    xi_m_diff_new[i + n] = xi_m_diff[i + 1 + n];
                    C_new.at(i) = C(i + 1);
                }
            }

            n--;
            xi_m_diff = xi_m_diff_new;
            C = C_new;

            logger.log("Removed at " + std::to_string(max_at) + " out of " + std::to_string(n) + ", " + std::to_string(max_diff));
        } else {
            logger.log("Max diff at " + std::to_string(max_at) + " out of " + std::to_string(n) + ", " + std::to_string(max_diff));
            break;
        }
    }

    logger.logv("Sat count after diff filtering (n)", n);

    // Сглаживание измерений
    for (size_t i = 0; i < n; ++i) {
        xi_m_diff[i] = 1 / Tp * xi_m_diff[i] + (1 - 1 / Tp) * xi_m_diff[i + n];
    }

    double trace;
    try {
        trace = sqrt((C0.transpose() * C0).inverse().trace());
        logger.logv("Trace", trace);
    } catch (const std::runtime_error&) {
        trace = C0_trace_threshold * 2;
        logger.log("Trace failed");
    }

    if (trace > C0_trace_threshold) {
        W = B1.transpose() * lambda * W * lambda * B1;
        dfs.x = x_est;
        dfs.dx = dx_est;
        
        ++consecutive_model_steps;
        since_last_model_step = 0;

        dfs_prev = dfs;
        solution.position = dfs.x;
        solution.is_solved = false;
        solution.failure_type = 'T';

        logger.log("Pure modeling");
        logger.log("Bad trace");
        logger.logv("Model steps", consecutive_model_steps);
        
        return solution;
    }

    if (W.norm() > W_norm_threshold) {
        W = zero(6, 6);
        // dfs.x = coarse;
        dfs.x = x_est;
        dfs.dx = dx_est;

        ++consecutive_model_steps;
        since_last_model_step = 0;

        dfs_prev = dfs;
        solution.position = dfs.x;
        solution.is_solved = false;
        solution.failure_type = 'W';

        logger.log("Pure modeling");
        logger.logv("Reseted filter", W.norm());
        logger.logv("Model steps", consecutive_model_steps);

        return solution;
    }

    std::vector<double> P = C0.transpose() * xi_m_diff;
    W = B1.transpose() * lambda * W * lambda * B1 + C0.transpose() * C0;

    std::vector<double> d_xi(6);
    try {
        d_xi = W.inverse() * P;
        is_model_step = false;
    } catch (const std::runtime_error& runtime_error) {
        is_model_step = true;
        logger.log("W is non-invertible: " + std::string(runtime_error.what()));
    }
    std::vector<double> xi = xi_est + d_xi;

    for (size_t i = 0; i < 3; ++i) {
        dfs.x[i] = xi[i];
        dfs.dx[i] = xi[i + 3];
    }

    logger.logv("W norm after", W.norm());
    logger.logv("Passive error", abs(dfs.pas_pos - pas.get_true_state_iterator(time)->position));
    logger.logv("dx error", abs(dfs.dx - (get_true_state_iterator(time)->position - get_true_state_iterator(time - dt)->position)));
    logger.logv("x error", abs(dfs.x - get_true_state_iterator(time)->position));

    if (is_model_step) {
        ++consecutive_model_steps;
        since_last_model_step = 0;
    } else {
        consecutive_model_steps = 0;
        ++since_last_model_step;
    }

    logger.logv("Is model step", is_model_step);
    logger.logv("Time since last model step", since_last_model_step);
    logger.logv("Consecutive model steps", consecutive_model_steps);

    if (since_last_model_step < model_steps_relaxation_threshold) {
        solution.is_solved = false;
        solution.failure_type = 'M';
        logger.log("Too soon since last model step");
    } else {
        solution.is_solved = true;
    }

    dfs_prev = dfs;
    solution.position = dfs.x;
    return solution;
}

std::vector<double> SatNavRel::estimate_dx() {
    double r = abs(dfs_prev.pas_pos);
    double omega2 = earth_mu / (r * r * r);
    Matrix D = 3 * tensor(dfs_prev.pas_pos, dfs_prev.pas_pos) / (r * r);

    return (E3 + 2 * Omega * dt) * dfs_prev.dx +
           -Omega * Omega * dfs_prev.x * dt * dt + 
           -omega2 * (E3 - D) * dfs_prev.x * dt * dt;
}

Matrix SatNavRel::calculate_B1() {
    double r = abs(dfs_prev.pas_pos);
    double omega2 = earth_mu / (r * r * r);
    Matrix D = 3 * tensor(dfs_prev.pas_pos, dfs_prev.pas_pos) / (r * r);
    Matrix derivative = -(Omega * Omega + omega2 * (E3 - D)) * dt * dt;

    Matrix B1(6, 6);
    B1.insert(E3, 0, 0);
    B1.insert(-E3, 0, 3);
    B1.insert(-derivative, 3, 0);
    B1.insert(E3 + derivative - 2 * Omega * dt, 3, 3);
    return B1;
}

double SatNavRel::calculate_delta(const std::vector<double>& L, const std::vector<double>& dX, const std::vector<double>& V) {
    double Ln = abs(L);
    double LdX = dot(L, dX);
    double dXdX = dot(dX, dX);
    double VdX = dot(V, dX);
    return dXdX / (2 * Ln) - LdX * LdX / (2 * Ln * Ln * Ln) - LdX * LdX * LdX / (2 * Ln * Ln * Ln * Ln * Ln) + LdX * dXdX / (2 * Ln * Ln * Ln) +
           VdX / c - VdX * VdX / (Ln * c * c) + VdX * VdX * VdX / (2 * Ln * Ln * c * c * c) - VdX * dXdX / (2 * Ln * Ln * c);
}

std::vector<double> SatNavRel::get_coarse_position(SatType sat_type, unsigned time) {
    SatNav* sat;

    switch (sat_type) {
    case SatType::PASSIVE:
        sat = &pas;
        break;
    case SatType::ACTIVE:
        sat = &act;
        break;
    }

    if (sat->get_solution_state_iterator(time)->is_solved) {
        logger.log("Coarse found");
        return sat->get_solution_state_iterator(time)->position;
    } else {
        logger.log("Coarse not found");
        return sat->get_true_state_iterator(time)->position;
    }
}

const std::vector<State>& SatNavRel::get_true_states() const {
    return true_states;
}

const std::vector<SolutionState>& SatNavRel::get_solution_states() const {
    return solution_states;
}

std::vector<State>::const_iterator SatNavRel::get_true_state_iterator(unsigned time) const {
    size_t n = true_states.size();

    std::ptrdiff_t low = 0;
    std::ptrdiff_t high = static_cast<std::ptrdiff_t>(n) - 1;
    while (low <= high) {
        std::ptrdiff_t mid = low + (high - low) / 2;

        if (true_states[static_cast<size_t>(mid)].time == time) {
            return true_states.cbegin() + mid;
        }

        if (true_states[static_cast<size_t>(mid)].time < time) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }

    return true_states.cend();
}