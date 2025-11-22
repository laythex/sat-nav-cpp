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

    void solve_relative(char et = '0', int ti = 0, int tf = 0);

    const State& get_true_state_at(int time) const;
    const std::vector<SolutionState>& get_solution_states() const;

    // должно быть приватным
    SatNav& pas;
    SatNav& act;

private:
    RefinedMeasurement refine_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act);
    SolutionState calculate_solution(const RefinedMeasurementGroupped& ref_mg) const;

    std::vector<State>::const_iterator get_true_state_iterator(int time) const;

    std::vector<State> true_states;
    std::vector<SolutionState> solution_states;
    std::vector<RefinedMeasurementGroupped> refined_measurements_groupped;

    char error_type;
};
