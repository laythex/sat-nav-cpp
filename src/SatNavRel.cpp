#include "SatNavRel.hpp"
#include "SatNav.hpp"

SatNavRel::SatNavRel(SatNav& passive, SatNav& active) : pas(passive), act(active) {}

void SatNavRel::solve_relative(char et = '0', int ti = 0, int tf = 0) {

    error_type = et;

    for (const auto& raw_mg_pas : pas.raw_measurements_groupped) {

        if ((ti > 0 || tf > 0) && (raw_mg_pas.time < ti || raw_mg_pas.time > tf)) continue;

        RawMeasurementGroupped raw_mg_act = act.get_raw_measurement_groupped_at(raw_mg_pas.time);

        if (raw_mg_act.time == -1) continue;

        std::vector<RefinedMeasurement> ref_ms(32);

        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            unsigned prn_index = prn_id - 1;

            if (!pas.check_raw(raw_mg_pas.raw_measurements[prn_index]) || !act.check_raw(raw_mg_act.raw_measurements[prn_index])) continue;

            RefinedMeasurement ref_m = refine_raw(raw_mg_pas.raw_measurements[prn_index], raw_mg_act.raw_measurements[prn_index]);

            ref_ms[prn_index] = ref_m;
        }

        RefinedMeasurementGroupped ref_mg = {raw_mg_pas.time, ref_ms};
        refined_measurements_groupped.push_back(ref_mg);

        SolutionState solution = calculate_solution(ref_mg);

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

        // solution_states.push_back(solution);
    }
   
    error_type = '0';
}

RefinedMeasurement SatNavRel::refine_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act) {
    double pseudorange_delta = raw_m_pas.L1_range - raw_m_act.L1_range;
    // Вычисляем положение НКА очень грубо
    // Может вычесть хотя бы ошибку часов НКА?
    std::vector<double> gps_position = pas.handler.get_state(raw_m_pas.prn_id, raw_m_pas.time).position;

    return {true, raw_m_pas.time, raw_m_pas.prn_id, pseudorange_delta, 0, gps_position};
}

SolutionState SatNavRel::calculate_solution(const RefinedMeasurementGroupped& ref_mg) const {
    SolutionState solution;
    solution.time = ref_mg.time;

    std::vector<double> dPR;
    std::vector<std::vector<double>> X;
    size_t n = 0;

    for (const auto& ref_m : ref_mg.refined_measurements) {
        if (ref_m.is_present) {
            dPR.push_back(ref_m.pseudorange);
            X.push_back(ref_m.gps_position);
            n++;
        }
    }

    std::vector<double> U(n);
    Matrix B(n, 3);

    for (size_t i = 0; i < n; i++) {
        size_t j = i < n - 1 ? i + 1 : 0;

        U[i] = 0;
        for (size_t k = 0; k < 3; k++) {
            B.at(i, k) = 0;
        }
    }
}
