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
    Matrix calculate_B1();
    double calculate_delta(const std::vector<double>& L, const std::vector<double>& dX, const std::vector<double>& V);
    std::vector<double> get_coarse(unsigned time);

    std::vector<State> true_states;
    std::vector<SolutionState> solution_states;
    std::vector<RefinedMeasurementGroupped> refined_measurements_groupped;

    Logger logger;

    char error_type;

    const double c = 2.99792458e8;
    const double earth_rotation_rate = 7.2921151467e-5;
    const double earth_mu = 3.986004418e14;

    const double CN0_min_threshold = 25;
    const double CN0_max_threshold = 65;
    const double C0_trace_threshold = 5;
    const double meas_diff_threshold = 10;
    const double W_norm_threshold = 500;
    const unsigned model_steps_meas_diff_threshold = 5;
    const unsigned model_steps_relaxation_threshold = 3;

    const bool is_pure_modeling_mode = false;
    const unsigned dt = 10;
    const double Tx = 10000;
    const double Tv = Tx / 1000;
    const double Tp = 2;
    Matrix lambda = {{{Tx / (Tx + 1), 0, 0, 0, 0, 0}, {0, Tx / (Tx + 1), 0, 0, 0, 0}, {0, 0, Tx / (Tx + 1), 0, 0, 0},
                      {0, 0, 0, Tv / (Tv + 1), 0, 0}, {0, 0, 0, 0, Tv / (Tv + 1), 0}, {0, 0, 0, 0, 0, Tv / (Tv + 1)}}};
    Matrix E3 = identity(3);
    Matrix Omega = {{{0, earth_rotation_rate, 0}, 
                     {-earth_rotation_rate, 0, 0}, 
                     {0, 0, 0}}};

    int df_state = 0;
    bool is_model_step;
    unsigned consecutive_model_steps = 0;
    unsigned since_last_model_step = 0;
    Matrix W = zero(6, 6);

    // Предыдущее состояние системы для динамической фильтрации
    DynamicFilterState dfs_prev;
};
