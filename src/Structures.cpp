#include "Structures.hpp"

State::State() : time(0), position(3, 0), velocity(3, 0) {}
State::State(int time, 
             const std::vector<double>& position,
             const std::vector<double>& velocity) :
             time(time), position(position), velocity(velocity) {}

GPSState::GPSState() : State(), relativistic_error(0) {}
GPSState::GPSState(int time, 
                   const std::vector<double>& position,
                   const std::vector<double>& velocity,
                   double relativistic_error) :
                   State(time, position, velocity), 
                   relativistic_error(relativistic_error) {}

SolutionState::SolutionState() : State(), is_solved(false), failure_type('0'), GDOP(0) {}
SolutionState::SolutionState(int time, 
                             const std::vector<double>& position,
                             const std::vector<double>& velocity,
                             bool is_solved, char failure_type, double GDOP) :
                             State(time, position, velocity),
                             is_solved(is_solved), failure_type(failure_type), GDOP(GDOP) {}

RawMeasurement::RawMeasurement() : is_present(false),
                                   time(0), prn_id(0), 
                                   L1_range(0), L2_range(0),
                                   L1_phase(0), L2_phase(0),
                                   L1_SNR(0), L2_SNR(0),
                                   qualflg(0) {}
RawMeasurement::RawMeasurement(bool is_present,
                               int time, unsigned prn_id,
                               double L1_range, double L2_range,
                               double L1_phase, double L2_phase,
                               unsigned short L1_SNR, unsigned short L2_SNR,
                               unsigned qualflg) : 
                               is_present(is_present),
                               time(time), prn_id(prn_id), 
                               L1_range(L1_range), L2_range(L2_range), 
                               L1_phase(L1_phase), L2_phase(L2_phase),
                               L1_SNR(L1_SNR), L2_SNR(L2_SNR),
                               qualflg(qualflg) {}

RefinedMeasurement::RefinedMeasurement() : is_present(false), 
                                           time(0), prn_id(0), pseudorange(0), carrier_phase(0),
                                           gps_position(3, 0), gps_velocity(3, 0) {}
RefinedMeasurement::RefinedMeasurement(bool is_present, 
                                       int time, unsigned prn_id, double pseudorange, double carrier_phase,
                                       const std::vector<double>& gps_position, const std::vector<double>& gps_velocity) : 
                                       is_present(is_present), 
                                       time(time), prn_id(prn_id), pseudorange(pseudorange), carrier_phase(carrier_phase),
                                       gps_position(gps_position), gps_velocity(gps_velocity) {}

RawMeasurementGroupped::RawMeasurementGroupped() : time(-1), raw_measurements(32) {}
RawMeasurementGroupped::RawMeasurementGroupped(int time,
                                               const std::vector<RawMeasurement>& raw_measurements) :
                                               time(time), raw_measurements(raw_measurements) {}


RefinedMeasurementGroupped::RefinedMeasurementGroupped() : time(-1), refined_measurements(32) {}
RefinedMeasurementGroupped::RefinedMeasurementGroupped(int time,
                                                       const std::vector<RefinedMeasurement>& refined_measurements) :
                                                       time(time), refined_measurements(refined_measurements) {}

DynamicFilterState::DynamicFilterState() : time(0), x(3), dx(3), mask(32, false), xi_m_1(32), xi_m_2(32), pas_pos(3), C_rows(32), is_first(false), is_second(false), is_failed(false) {}
DynamicFilterState::DynamicFilterState(int time, 
                                       const std::vector<double>& x, const std::vector<double>& dx,
                                       const std::vector<bool>& mask, const std::vector<double>& xi_m_1, const std::vector<double>& xi_m_2,
                                       const std::vector<double>& pas_pos, const std::vector<std::vector<double>>& C_rows,
                                       bool is_first, bool is_second, bool is_failed) : 
                                       time(time), x(x), dx(dx), 
                                       mask(mask), xi_m_1(xi_m_1), xi_m_2(xi_m_2),
                                       pas_pos(pas_pos), C_rows(C_rows),
                                       is_first(is_first), is_second(is_second), is_failed(is_failed) {}
