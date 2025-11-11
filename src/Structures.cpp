#include "Structures.hpp"

State::State() : time(0), position(3, 0), velocity(3, 0) {}

GPSState::GPSState() : State(), relativistic_error(0) {}

SolutionState::SolutionState() : State(), is_solved(false), failure_type('0'), GDOP(0) {}

RawMeasurement::RawMeasurement() : time(0), prn_id(0), 
                               L1_range(0), L2_range(0),
                               L1_phase(0), L2_phase(0),
                               L1_SNR(0), L2_SNR(0),
                               qualflg(0) {}

RefinedMeasurement::RefinedMeasurement() : time(0), prn_id(0), pseudorange(0), gps_position(3, 0) {}

RawMeasurementGroupped::RawMeasurementGroupped() : time(0), raw_measurements(32) {}

RefinedMeasurementGroupped::RefinedMeasurementGroupped() : time(0), refined_measurements(32) {}
