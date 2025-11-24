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

    void solve_relative(char et = '0', double ti = 0, double tf = 0);

    const std::vector<State>& get_true_states() const;
    const std::vector<SolutionState>& get_solution_states() const;

    std::vector<State>::const_iterator get_true_state_iterator(int time) const;

    // должно быть приватным?
    SatNav& pas;
    SatNav& act;

private:
    RefinedMeasurement refine_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act);
    SolutionState calculate_solution(const RefinedMeasurementGroupped& ref_mg, const std::vector<double>& approximation);
    std::vector<SolutionState>::const_reverse_iterator get_last_solved() const;

    std::vector<State> true_states;
    std::vector<SolutionState> solution_states;
    std::vector<RefinedMeasurementGroupped> refined_measurements_groupped;

    Logger logger;

    char error_type;

    const double c = 2.99792458e8;

    const double GDOP0 = 5;
    const double approximation_threshold = 10;
};
