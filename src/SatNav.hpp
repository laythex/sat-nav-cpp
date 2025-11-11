#pragma once

#include <iostream>

#include <vector>

#include "LinAlg.hpp"
#include "Structures.hpp"
#include "GPSHandler.hpp"
#include "DataParser.hpp"


class SatNav {

public:
    SatNav(const std::string& gnv_filename, const std::string& gps_filename, const GPSHandler& handler);
    void solve(unsigned ti = 0, unsigned tf = 0, char error_type = '0');

private:
    RefinedMeasurement refine_raw(const RawMeasurement& raw_m);
    RefinedMeasurement apply_clock_and_relativistic_errors(const RawMeasurement& raw_m, unsigned frequency);
    SolutionState calculate_solution(const RefinedMeasurementGroupped& ref_mg) const;

    // std::vector<double> get_true_position(unsigned user_time);
    // std::vector<double> get_true_velocity(unsigned user_time);
    
    GPSHandler handler;

    std::vector<State> true_states;
    std::vector<SolutionState> solution_states;
    std::vector<RawMeasurementGroupped> raw_measurements_groupped;
    std::vector<RefinedMeasurementGroupped> refined_measurements_groupped;

    double c = 2.99792458e8;
    double earth_rotation_rate = 7.2921151467e-5;
    double gamma = 1.64694444444;
    double C1 = 2.54572778016;
    double C2 = -1.54572778016;

    double GDOP0 = 5;
    double mask_angle = 5;
    unsigned fadeout_time = 20;
    double SNR_threshold = 30;
    double hatch_constant = 0.1;

    char error_type = '0';
};
