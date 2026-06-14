#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "Config.hpp"
#include "Constants.hpp"
#include "DataParser.hpp"
#include "GPSHandler.hpp"
#include "Logger.hpp"
#include "SatNavRel.hpp"
#include "SatNavUtils.hpp"
#include "Structures.hpp"

class SatNavRel;

class SatNav {

public:
    SatNav(const Date& date, const GRACE sat, const GPSHandler& handler);
    SatNav(const Date& date, const SWARM sat, const GPSHandler& handler);
    SatNav(const SatNav& sn);

    void solve(IntendedError error = IntendedError::NONE);

    const std::vector<State>& get_true_states() const;
    const std::vector<SolutionState>& get_solution_states() const;
    const std::vector<RawMeasurementGroupped>& get_raw_measurements_groupped() const;
    const std::vector<RefinedMeasurementGroupped>& get_refined_measurements_groupped() const;
    const std::vector<AccelerationMeasurement>& get_acceleration_measurements() const;

    std::vector<State>::const_iterator get_true_state_iterator(unsigned time) const;
    std::vector<SolutionState>::const_iterator get_solution_state_iterator(unsigned time) const;
    std::vector<RawMeasurementGroupped>::const_iterator get_raw_measurement_groupped_iterator(unsigned time) const;

    Date get_date() const;

    friend class SatNavRel;

private:
    SolutionState find_solution(const RawMeasurementGroupped& raw_mg);

    std::vector<size_t> check_raw_groupped(const RawMeasurementGroupped& raw_mg);
    MeasurementStatus check_raw(const RawMeasurement& raw_m) const;

    RefinedMeasurementGroupped refine_raw_groupped(const RawMeasurementGroupped& raw_mg,
                                                   const std::vector<size_t>& present_prns) const;
    RefinedMeasurement refine_raw(const RawMeasurement& raw_m) const;

    bool check_fadein(const RawMeasurement& raw_m) const;
    RefinedMeasurement hatch_filter(const RefinedMeasurement& raw_m) const;

    Date date;
    GPSHandler handler;
    Logger logger;

    std::vector<State> true_states;
    std::vector<RawMeasurementGroupped> raw_measurements_groupped;
    std::vector<AccelerationMeasurement> acceleration_measurements;
    std::vector<SolutionState> solution_states;
    std::vector<RefinedMeasurementGroupped> refined_measurements_groupped;

    IntendedError error = IntendedError::NONE;
};
