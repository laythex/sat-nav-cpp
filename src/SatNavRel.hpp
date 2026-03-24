#pragma once

#include <iostream>

#include <vector>

#include "LinAlg.hpp"
#include "Structures.hpp"
#include "GPSHandler.hpp"
#include "DataParser.hpp"
#include "Logger.hpp"
#include "SatNav.hpp"
#include "Propagator.hpp"

class SatNav;

class SatNavRel {

public:
    SatNavRel(SatNav& passive, SatNav& active);
    SatNavRel(const SatNavRel& sn);

    void solve_relative(char et = '0', int ti = 0, int tf = 0);

    const std::vector<State>& get_true_states() const;
    const std::vector<SolutionState>& get_solution_states() const;

    std::vector<State>::const_iterator get_true_state_iterator(int time) const;

    // должно быть приватным?
    SatNav& pas;
    SatNav& act;

private:
    // Проверяем что с сырыми измерениями все ок
    bool check_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act);

    // Обрабатываем сырые измерения
    RefinedMeasurement refine_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act);

    // Находим вектор состояния
    SolutionState calculate_solution(const RefinedMeasurementGroupped& ref_mg);

    std::vector<double> estimate_dx();
    std::vector<double> estimate_dpp();
    Matrix calculate_B1(const std::vector<double>& pas_pos);

    std::vector<State> true_states;
    std::vector<SolutionState> solution_states;
    std::vector<RefinedMeasurementGroupped> refined_measurements_groupped;

    Logger logger;

    char error_type;

    const double c = 2.99792458e8;
    const double earth_rotation_rate = 7.2921151467e-5;
    const double earth_mu = 3.986004418e14;
    const double CN0_constant = 3.01029995664;

    const double GDOP0 = 3;
    const double CN0_min = 28;
    const double CN0_max = 65;

    // Глобальные перменные динамического фильтра
    bool df_started = false;
    bool df_half_started = false;
    double T_x = 100;
    double T_v = 100;
    double T_p = 1;
    Matrix lambda = Matrix(6, 6, 0.0);
    Matrix E3 = identity(3);
    Matrix Omega = {{{0, earth_rotation_rate, 0}, {-earth_rotation_rate, 0, 0}, {0, 0, 0}}};
    Matrix W = Matrix(6, 6, 0.0);

    // Предыдущее состояние системы для динамической фильтрации
    DynamicFilterState dfs_prev;
};
