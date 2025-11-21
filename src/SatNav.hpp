#pragma once

#include <iostream>

#include <vector>

#include "LinAlg.hpp"
#include "Structures.hpp"
#include "GPSHandler.hpp"
#include "DataParser.hpp"
#include "Logger.hpp"
#include "SatNavRel.hpp"

class SatNavRel;

class SatNav {

public:
    SatNav(const std::string& gnv_filename, const std::string& gps_filename, const GPSHandler& handler);
    SatNav(const SatNav& sn);

    void solve(char et = '0', int ti = 0, int tf = 0);

    const State& get_true_state_at(int time) const;
    const SolutionState& get_solution_state_at(int time) const;
    const RawMeasurementGroupped& get_raw_measurement_groupped_at(int time) const;

    const std::vector<State>& get_true_states() const;
    const std::vector<SolutionState>& get_solution_states() const;
    const std::vector<RawMeasurementGroupped>& get_raw_measurements_groupped() const;
    const std::vector<RefinedMeasurementGroupped>& get_refined_measurements_groupped() const;

    friend class SatNavRel;

private:
    bool check_raw(const RawMeasurement& raw_m);
    RefinedMeasurement refine_raw(const RawMeasurement& raw_m);
    RefinedMeasurement apply_clock_and_relativistic_errors(const RawMeasurement& raw_m, unsigned frequency);
    SolutionState calculate_solution(const RefinedMeasurementGroupped& ref_mg) const;
    bool check_fading(const RawMeasurement& raw_m);
    std::vector<unsigned> check_low(const SolutionState& solution, const RefinedMeasurementGroupped& ref_mg);
    RefinedMeasurement hatch_filter(const RefinedMeasurement& raw_m);

    std::vector<State>::const_iterator get_true_state_iterator(int time) const;
    std::vector<SolutionState>::const_iterator get_solution_state_iterator(int time) const;
    std::vector<RawMeasurementGroupped>::const_iterator get_raw_measurement_groupped_iterator(int time) const;

    GPSHandler handler;

    std::vector<State> true_states;
    std::vector<SolutionState> solution_states;
    std::vector<RawMeasurementGroupped> raw_measurements_groupped;
    std::vector<RefinedMeasurementGroupped> refined_measurements_groupped;

    Logger logger;

    double c = 2.99792458e8;
    double earth_rotation_rate = 7.2921151467e-5;
    double C1 = 2.54572778016;
    double C2 = -1.54572778016;

    double GDOP0 = 5;
    double mask_angle = 5;
    int fadeout_time = 20;
    double SNR_threshold = 30;
    double hatch_constant = 1.0 / 30.0;

    char error_type = '0';
};
