#include "SatNavRel.hpp"
#include "SatNav.hpp"

SatNavRel::SatNavRel(SatNav& passive, SatNav& active) : pas(passive), act(active) {
    pas.solve(); // решать одновременно?
    // act.solve();

    for (const State& ts_pas : pas.get_true_states()) {

        auto ts_act_it = act.get_true_state_iterator(ts_pas.time);
        if (ts_act_it == act.get_true_states().end()) continue;
        State ts_act = *ts_act_it;

        true_states.push_back({ts_pas.time, ts_act.position - ts_pas.position, ts_act.velocity - ts_pas.velocity});
    }
}

SatNavRel::SatNavRel(const SatNavRel& sn) : pas(sn.pas), act(sn.act),
                                            true_states(sn.true_states) {}

void SatNavRel::solve_relative(char et, double ti, double tf) {

    logger.log("Beginning to solve...");

    error_type = et;

    size_t raw_mg_cnt = 0;
    size_t raw_mg_size = pas.raw_measurements_groupped.size();
    for (const auto& raw_mg_pas : pas.raw_measurements_groupped) {

        if (ti > 0 || tf > 0) {
            double t = static_cast<double>(raw_mg_cnt) / static_cast<double>(raw_mg_size);
            if (t < ti || t >= tf) continue;
            raw_mg_cnt++;
        }

        auto raw_mg_act_it = act.get_raw_measurement_groupped_iterator(raw_mg_pas.time);
        if (raw_mg_act_it == act.raw_measurements_groupped.end()) continue;
        RawMeasurementGroupped raw_mg_act = *raw_mg_act_it;

        if (raw_mg_act.time == -1) continue;

        logger.log("Current time = " + std::to_string(raw_mg_pas.time));

        std::vector<RefinedMeasurement> ref_ms(32);

        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            unsigned prn_index = prn_id - 1;

            if (!pas.check_raw(raw_mg_pas.raw_measurements[prn_index]) || !act.check_raw(raw_mg_act.raw_measurements[prn_index])) continue;

            logger.log("Refining " + std::to_string(prn_id));

            RefinedMeasurement ref_m = refine_raw(raw_mg_pas.raw_measurements[prn_index], raw_mg_act.raw_measurements[prn_index]);
            ref_ms[prn_index] = ref_m;
        }

        RefinedMeasurementGroupped ref_mg = {raw_mg_pas.time, ref_ms};
        refined_measurements_groupped.push_back(ref_mg);

        bool last_solved_found = false;
        std::vector<double> approximation = {0.0, 0.0, 0.0};

        auto ss_last_it = get_last_solved();
        if (ss_last_it != solution_states.rend()) {
            last_solved_found = true;
            approximation = ss_last_it->position;
        }

        SolutionState solution = calculate_solution(ref_mg, approximation);
        if (!last_solved_found || solution.time - ss_last_it->time > approximation_threshold) {
            solution = calculate_solution(ref_mg, solution.position);
        }

        logger.log("Is solved: " + std::to_string(solution.is_solved) + ' ' + solution.failure_type);

        if (solution.is_solved) {
            logger.log("GDOP: " + std::to_string(solution.GDOP));
            logger.log("Solution norm: " + std::to_string(abs(solution.position)));

            State ts = *get_true_state_iterator(raw_mg_pas.time);
            double error = abs(solution.position - ts.position);
            logger.log("Error: " + std::to_string(error));
            if (error > 5) logger.log("HIGH ERROR");
            // 187512340
            SolutionState ss_pas = *pas.get_solution_state_iterator(raw_mg_pas.time);
            State ts_pas = *pas.get_true_state_iterator(raw_mg_pas.time);
            double error_pas = abs(ss_pas.position - ts_pas.position);
            logger.log("Passive error: " + std::to_string(error_pas));
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
    double pseudorange_delta = raw_m_pas.L1_range - raw_m_act.L1_range;

    double clock_error = pas.handler.get_clock_error(raw_m_pas.prn_id, raw_m_pas.time);
    double delay = raw_m_pas.L1_range / 3e8 + clock_error;

    GPSState gs = pas.handler.get_state(raw_m_pas.prn_id, raw_m_pas.time - delay);
    
    return {true, raw_m_pas.time, raw_m_pas.prn_id, pseudorange_delta, 0, gs.position, gs.velocity};
}

SolutionState SatNavRel::calculate_solution(const RefinedMeasurementGroupped& ref_mg, const std::vector<double>& approximation) {
    SolutionState solution;
    solution.time = ref_mg.time;

    auto ss_pas_it = pas.get_solution_state_iterator(ref_mg.time);
    if (ss_pas_it == pas.get_solution_states().end()) {
        solution.is_solved = false;
        solution.failure_type = 'p';
        return solution;
    }

    SolutionState ss_pas = *ss_pas_it;
    
    if (!ss_pas.is_solved) {
        solution.is_solved = false;
        solution.failure_type = 'P';
        return solution;
    }

    std::vector<double> pr_pas_minus_pr_act;
    std::vector<std::vector<double>> gps_pos_minus_pas_pos;
    std::vector<double> deltas;

    size_t n = 0;

    for (const auto& ref_m : ref_mg.refined_measurements) {
        if (ref_m.is_present) {
            std::vector<double> X_gps_m_X_pas = ref_m.gps_position - ss_pas.position;

            pr_pas_minus_pr_act.push_back(ref_m.pseudorange);
            gps_pos_minus_pas_pos.push_back(X_gps_m_X_pas);

            double delta = ref_m.gps_velocity * approximation / (c * c) + 
                           (X_gps_m_X_pas * approximation) * (X_gps_m_X_pas * approximation) / (abs(X_gps_m_X_pas) * abs(X_gps_m_X_pas) * abs(X_gps_m_X_pas)) / 2 -
                           approximation * approximation / abs(X_gps_m_X_pas) / 2;

            deltas.push_back(delta);
            n++;
        }
    }

    std::vector<double> U(n);
    Matrix B(n, 3);

    for (size_t i = 0; i < n; i++) {
        size_t j = i < n - 1 ? i + 1 : 0;

        U[i] = pr_pas_minus_pr_act[i] - pr_pas_minus_pr_act[j] - (deltas[i] - deltas[j]);
        B.at(i) = gps_pos_minus_pas_pos[i] / abs(gps_pos_minus_pas_pos[i]) - 
                  gps_pos_minus_pas_pos[j] / abs(gps_pos_minus_pas_pos[j]);

    }

    Matrix B1(3, 3);
    try {
        B1 = (B.transpose() * B).inverse();
    } catch (std::runtime_error re) {
        solution.is_solved = false;
        solution.failure_type = 'I';
        return solution;
    }

    double GDOP = sqrt(B1.trace());
    solution.GDOP = GDOP;

    if (GDOP > GDOP0) {
        solution.is_solved = false;
        solution.failure_type = 'G';
        return solution;
    }

    std::vector<double> dX = B1 * B.transpose() * U;

    solution.position = dX;
    solution.is_solved = true;

    return solution;
}

std::vector<SolutionState>::const_reverse_iterator SatNavRel::get_last_solved() const {
    for (auto ss_it = solution_states.rbegin(); ss_it != solution_states.rend(); ss_it++) {
        if (ss_it->is_solved) {
            return ss_it;
        }
    }

    return solution_states.rend();
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
