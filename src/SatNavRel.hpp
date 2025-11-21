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

private:
    RefinedMeasurement refine_raw(const RawMeasurement& raw_m_pas, const RawMeasurement& raw_m_act);
    SolutionState calculate_solution(const RefinedMeasurementGroupped& ref_mg) const;

    SatNav& pas;
    SatNav& act;

    std::vector<SolutionState> solution_states;
    std::vector<RefinedMeasurementGroupped> refined_measurements_groupped;

    char error_type;
};
