#pragma once

#include <iostream>
#include <vector>
#include <cstddef>

#include "LinAlg.hpp"
#include "Structures.hpp"
#include "GPSHandler.hpp"
#include "DataParser.hpp"
#include "Logger.hpp"
#include "SatNav.hpp"

class SatNav;

class SatNavRel {

public:
    SatNavRel(SatNav& passive, SatNav& active);
    SatNavRel(const SatNavRel& sn);

    void solve_relative(char et = '0', unsigned ti = 0, unsigned tf = 0);

    const std::vector<State>& get_true_states() const;
    const std::vector<SolutionState>& get_solution_states() const;

    std::vector<State>::const_iterator get_true_state_iterator(unsigned time) const;

    // должно быть приватным?
    SatNav& pas;
    SatNav& act;

private:
    void solve_separately(char et = '0', unsigned ti = 0, unsigned tf = 0);

    // Проверяем что с сырыми измерениями все ок
    bool check_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act);

    // Обрабатываем сырые измерения
    RefinedMeasurement refine_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act);

    // Находим вектор состояния
    SolutionState calculate_solution(const RefinedMeasurementGroupped& ref_mg);

    std::vector<double> estimate_dx();
    std::vector<double> estimate_dpp();
    Matrix calculate_B1();

    std::vector<State> true_states;
    std::vector<SolutionState> solution_states;
    std::vector<RefinedMeasurementGroupped> refined_measurements_groupped;

    Logger logger;

    char error_type;

    const double c = 2.99792458e8;
    const double earth_rotation_rate = 7.2921151467e-5;
    const double earth_mu = 3.986004418e14;
    const double CN0_constant = 3.01029995664;

    const double CN0_min_threshold = 25;
    const double CN0_max_threshold = 65;
    const double C0_trace_threshold = 8;
    const double xi_m_diff_threshold = 10;
    const double W_norm_threshold = 100;
    const unsigned model_steps_threshold = 10;
    const unsigned model_steps_relaxation = 3;

    const double T_x = 10000 * 0;
    const double T_v = T_x;
    const double T_p = 2;
    Matrix lambda = {{{T_x / (T_x + 1), 0, 0, 0, 0, 0}, {0, T_x / (T_x + 1), 0, 0, 0, 0}, {0, 0, T_x / (T_x + 1), 0, 0, 0},
                      {0, 0, 0, T_v / (T_v + 1), 0, 0}, {0, 0, 0, 0, T_v / (T_v + 1), 0}, {0, 0, 0, 0, 0, T_v / (T_v + 1)}}};
    Matrix E3 = identity(3);
    Matrix Omega = {{{0, earth_rotation_rate, 0}, {-earth_rotation_rate, 0, 0}, {0, 0, 0}}};

    int df_state = 0;
    unsigned model_steps = 0;
    unsigned since_model_steps = 0;
    Matrix W = zero(6, 6);

    // Предыдущее состояние системы для динамической фильтрации
    DynamicFilterState dfs_prev;
};
