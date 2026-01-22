#pragma once

#include <iostream>

#include <vector>

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

    void solve_relative(char et = '0', int ti = 0, int tf = 0);

    const std::vector<State>& get_true_states() const;
    const std::vector<SolutionState>& get_solution_states() const;

    std::vector<State>::const_iterator get_true_state_iterator(int time) const;
    std::vector<SolutionState>::const_iterator get_coarse_solution_state_iterator(int time) const;

    // должно быть приватным?
    SatNav& pas;
    SatNav& act;

private:
    // Обрабатываем сырые измерения (учитываем различные эффекты + считаем положение НКА)
    RefinedMeasurement refine_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act);

    // Находим вектор состояния
    SolutionState calculate_solution(const RefinedMeasurementGroupped& ref_mg);

    // Находим приблизительный вектор состояния
    State estimate_state(const SolutionState& ss_pas, const SolutionState& ss_last, int time);

    // Находим приблизительный вектор измерений

    std::vector<State> true_states;
    std::vector<SolutionState> solution_states;
    std::vector<SolutionState> coarse_solution_states;
    std::vector<RefinedMeasurementGroupped> refined_measurements_groupped;

    Logger logger;

    char error_type;

    const double c = 2.99792458e8;
    const double earth_rotation_rate = 7.2921151467e-5;

    const double GDOP0 = 5;

    // Постоянные динамического фильтра
    double T_x = 1;
    double T_v = 1;
    Matrix lambda = Matrix(6, 6, 0.0);
    Matrix W = Matrix(6, 6, 0.0);
};
