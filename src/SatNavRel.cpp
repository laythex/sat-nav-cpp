#include "SatNavRel.hpp"
#include "SatNav.hpp"

SatNavRel::SatNavRel(SatNav& passive, SatNav& active) : pas(passive), act(active) {}

SatNavRel::SatNavRel(const SatNavRel& sn) : pas(sn.pas), act(sn.act), true_states(sn.true_states) {}
                                        
void SatNavRel::solve_separately(char et, unsigned ti, unsigned tf) {
    pas.solve(et, ti, tf);
    act.solve(et, ti, tf);

    unsigned t0 = pas.get_true_states()[0].time;
    for (const State& ts_pas : pas.get_true_states()) {

        unsigned time = ts_pas.time;
        if (ti > 0 || tf > 0) {
            if ((time - t0) < ti || (time - t0) >= tf) continue;
        }

        auto ts_act_it = act.get_true_state_iterator(time);
        if (ts_act_it == act.get_true_states().end()) continue;
        State ts_act = *ts_act_it;

        true_states.push_back({time, ts_act.position - ts_pas.position, ts_act.velocity - ts_pas.velocity});
    }
}

void SatNavRel::solve_relative(char et, unsigned ti, unsigned tf) {

    solve_separately(et, ti, tf);

    logger.log("Starting to solve...");

    error_type = et;
    double err_avg = 0;

    unsigned t0 = pas.raw_measurements_groupped[0].time;
    for (const auto& raw_mg_pas : pas.raw_measurements_groupped) {

        unsigned time = raw_mg_pas.time;
        if (ti > 0 || tf > 0) {
            if ((time - t0) < ti || (time - t0) >= tf) continue;
        }

        auto raw_mg_act_it = act.get_raw_measurement_groupped_iterator(time);
        if (raw_mg_act_it == act.raw_measurements_groupped.end()) continue;
        RawMeasurementGroupped raw_mg_act = *raw_mg_act_it;

        logger.logv("Time", time);
        
        std::cout << "Time: " << time << std::endl;

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
            RefinedMeasurement ref_m = refine_raw(raw_m_pas, raw_m_act);
            ref_ms[prn_index] = ref_m;
        }

        // Группируем обработанные измерения со всех спутников
        RefinedMeasurementGroupped ref_mg = {time, ref_ms};
        refined_measurements_groupped.push_back(ref_mg);

        // Ищем решение
        SolutionState solution;
        solution = calculate_solution(ref_mg);
        solution_states.push_back(solution);

        logger.log("Is solved: " + std::to_string(solution.is_solved) + ' ' + solution.failure_type);
        auto ts = get_true_state_iterator(raw_mg_pas.time);
        double error = abs(solution.position - ts->position);
        logger.logv("Error norm", error);
        if (solution.is_solved) {
            err_avg += error;
            if (error > 3) {
                logger.log("HIGH ERROR");
            }
        }
        logger.lnbr();
    }
   
    error_type = '0';

    err_avg /= static_cast<double>(solution_states.size());
    logger.logv("Average error", err_avg);
    std::cout << "Avg error: " << err_avg << std::endl;

    logger.log("Finished solving");
}

bool SatNavRel::check_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act) {
    if (!raw_m_pas.is_present || !raw_m_act.is_present) {
        return false;
    }

    if (raw_m_pas.qualflg != 0 || raw_m_act.qualflg != 0) {
        return false;
    }

    if (error_type != 'S') {
        double L1_CN0_pas = 20 * log10(raw_m_pas.L1_SNR) - CN0_constant;
        double L2_CN0_pas = 20 * log10(raw_m_pas.L2_SNR) - CN0_constant;
        double L1_CN0_act = 20 * log10(raw_m_act.L1_SNR) - CN0_constant;
        double L2_CN0_act = 20 * log10(raw_m_act.L2_SNR) - CN0_constant;
        double min_CN0 = std::min(std::min(L1_CN0_pas, L2_CN0_pas), std::min(L1_CN0_act, L2_CN0_act));
        double max_CN0 = std::max(std::max(L1_CN0_pas, L2_CN0_pas), std::max(L1_CN0_act, L2_CN0_act));

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

    double clock_error = pas.handler_.get_clock_error(raw_m_pas.prn_id, raw_m_pas.time); // Лишнее?
    double delay = raw_m_pas.L1_range / c + clock_error;

    GPSState gs = pas.handler_.get_state(raw_m_pas.prn_id, raw_m_pas.time - delay);

    return {true, raw_m_pas.time, raw_m_pas.prn_id, pseudorange_delta, carrier_phase_delta, gs.position, gs.velocity};
}

SolutionState SatNavRel::calculate_solution(const RefinedMeasurementGroupped& ref_mg) {
    unsigned time = ref_mg.time;

    SolutionState solution;
    solution.time = time;

    DynamicFilterState dfs;
    dfs.time = time;

    // Находим решения пассивного и активного КА
    auto ss_pas_it = pas.get_solution_state_iterator(time);
    auto ss_act_it = act.get_solution_state_iterator(time);

    bool has_passive = ss_pas_it != pas.get_solution_states().end() && ss_pas_it->is_solved;
    bool has_active = ss_act_it != act.get_solution_states().end() && ss_act_it->is_solved;
    bool good_coarse_data = has_passive && has_active;

    if (df_state < 2 && !good_coarse_data) {
        df_state = 0;

        solution.is_solved = false;
        solution.failure_type = '-';
        return solution;
    }

    std::vector<double> pas_pos, act_pos;
    if (has_passive) {
        pas_pos = ss_pas_it->position;
    } else {
        pas_pos =  dfs_prev.pas_pos + estimate_dpp();
    }
    if (has_active) {
        act_pos = ss_act_it->position;
    } else {
        act_pos = pas_pos;
    }

    // продолжить тут

    good_coarse_data = true;
    pas_pos = pas.get_true_state_iterator(time)->position;
    act_pos = act.get_true_state_iterator(time)->position;
  
    std::vector<double> coarse = act_pos - pas_pos;
    dfs.pas_pos = pas_pos;

    // Проходим по всем присутствующим НКА
    for (const auto& ref_m : ref_mg.refined_measurements) {
        if (!ref_m.is_present) continue;

        std::vector<double> XiX = ref_m.gps_position - dfs.pas_pos;
        double Di = abs(XiX);
        double delta = 0.5 * (XiX * coarse) * (XiX * coarse) / (Di * Di * Di) +
                       -0.5 * coarse * coarse / Di;

        size_t prn_index = ref_m.prn_id - 1;
        dfs.mask[prn_index] = true;
        dfs.xi_m_1[prn_index] = ref_m.pseudorange - delta;
        dfs.xi_m_2[prn_index] = ref_m.carrier_phase;
        dfs.C_rows[prn_index] = XiX / Di;
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

    // Строим оценку вектора состояния путем прогноза
    std::vector<double> dx_est = estimate_dx();
    std::vector<double> x_est = dfs_prev.x + dx_est;
    dfs.x_est = x_est; // Запоминаем для следующей итерации

    // Строим матрицу B1
    Matrix B1 = calculate_B1();

    if (!good_coarse_data) {
        W = B1.transpose() * lambda * W * lambda * B1;
        dfs.x = x_est;
        dfs.dx = dx_est;
        
        ++model_steps;
        since_model_steps = 0;

        dfs_prev = dfs;
        solution.position = dfs.x;
        solution.is_solved = false;
        solution.failure_type = 'C';

        logger.log("Pure modeling");
        logger.log("Bad coarse data");
        logger.logv("Model steps", model_steps);
        
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

    // Строим оценку вектора измерений (3.55)
    std::vector<double> x_m_est = C * x_est;
    std::vector<double> dx_m_est = dC * dfs_prev.x_est + C * dx_est;

    logger.logv("Passive error", abs(dfs.pas_pos - pas.get_true_state_iterator(time)->position));
    logger.logv("Previous error", abs(dfs_prev.x - get_true_state_iterator(time - 10)->position));
    logger.logv("Model error", abs(x_est - get_true_state_iterator(time)->position));

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

    // 3.68
    while (n > 0) {
        if (model_steps > model_steps_threshold) {
            logger.logv("Model steps exceed threshold, no filtering", model_steps);
            break;
        }

        size_t max_at = find_max_abs(xi_m_diff);
        double max_diff = std::abs(xi_m_diff[max_at]);

        if (max_diff > xi_m_diff_threshold) {
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
        xi_m_diff[i] = 1 / T_p * xi_m_diff[i] + (1 - 1 / T_p) * xi_m_diff[i + n];
    }

    Matrix C0 = zero(n * 2, 6);
    C0.insert(C, 0, 0);
    C0.insert(C, n, 3);

    Matrix C0_Gram = C0.transpose() * C0;
    bool good_trace = false;
    try {
        double trace = sqrt(C0_Gram.inverse().trace());
        good_trace = trace < C0_trace_threshold;
        logger.logv("Trace", trace);
    } catch (const std::runtime_error&) {
        logger.log("Trace failed");
    }

    if (!good_trace) {
        W = B1.transpose() * lambda * W * lambda * B1;
        dfs.x = x_est;
        dfs.dx = dx_est;
        
        ++model_steps;
        since_model_steps = 0;

        dfs_prev = dfs;
        solution.position = dfs.x;
        solution.is_solved = false;
        solution.failure_type = 'T';

        logger.log("Pure modeling");
        logger.log("Bad trace");
        logger.logv("Model steps", model_steps);
        
        return solution;
    }

    if (W.norm() > W_norm_threshold) {
        W = zero(6, 6);
        dfs.x = x_est;
        dfs.dx = dx_est;

        ++model_steps;
        since_model_steps = 0;

        dfs_prev = dfs;
        solution.position = dfs.x;
        solution.is_solved = false;
        solution.failure_type = 'W';

        logger.log("Pure modeling");
        logger.logv("Reseted W", W.norm());
        logger.logv("Model steps", model_steps);

        return solution;
    }

    // Динамическая фильтрация
    std::vector<double> P = C0.transpose() * xi_m_diff;
    W = B1.transpose() * lambda * W * lambda * B1 + C0_Gram;

    logger.logv("W norm after", W.norm());

    // Проверка обратимости W
    std::vector<double> d_xi(6);
    try {
        d_xi = W.inverse() * P;
    } catch (const std::runtime_error& runtime_error) {
        logger.log("W is non-invertible: " + std::string(runtime_error.what()));
    }
    std::vector<double> xi = xi_est + d_xi;

    // Сохраняем на следующую итерацию
    for (size_t i = 0; i < 3; ++i) {
        dfs.x[i] = xi[i];
        dfs.dx[i] = xi[i + 3];
    }

    model_steps = 0;
    ++since_model_steps;

    dfs_prev = dfs;
    solution.position = dfs.x;

    logger.logv("Time since last model step", since_model_steps);
    if (since_model_steps < model_steps_relaxation) {
        solution.is_solved = false;
        solution.failure_type = 'M';
        logger.log("Too soon since last model step");
    } else {
        solution.is_solved = true;
    }

    return solution;
}

std::vector<double> SatNavRel::estimate_dpp() {
    double dt = 10;
    double r = abs(dfs_prev.pas_pos);
    double omega2 = earth_mu / (r * r * r);

    return (E3 + 2 * Omega * dt) * dfs_prev.d_pas_pos - 
           Omega * Omega * dfs_prev.pas_pos * dt * dt - 
           omega2 * dfs_prev.pas_pos * dt * dt;
}

std::vector<double> SatNavRel::estimate_dx() {
    double dt = 10;
    double r = abs(dfs_prev.pas_pos);
    double omega2 = earth_mu / (r * r * r);
    Matrix D = 3 * tensor(dfs_prev.pas_pos, dfs_prev.pas_pos) / (r * r);

    return (E3 + 2 * Omega * dt) * dfs_prev.dx +
           -Omega * Omega * dfs_prev.x * dt * dt + 
           -omega2 * (E3 - D) * dfs_prev.x * dt * dt;
}

Matrix SatNavRel::calculate_B1() {
    double dt = 10;
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