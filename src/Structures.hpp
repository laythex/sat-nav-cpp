#pragma once

#include <vector>

struct Ephemeris {
    unsigned prn_id;
    int t_oe;

    double a_f0, a_f1, a_f2;
    double M_0, delta_n, e, A_sqrt, Omega_0, i_0, omega, Omega_dot, IDOT;
    double C_uc, C_us, C_rc, C_rs, C_ic, C_is;
};

struct State {
    int time;
    std::vector<double> position;
    std::vector<double> velocity;

    State();
    State(int time, 
          const std::vector<double>& position,
          const std::vector<double>& velocity);
};

struct GPSState : State {
    double relativistic_error;

    GPSState();
    GPSState(int time, 
             const std::vector<double>& position,
             const std::vector<double>& velocity,
             double relativistic_error);
};

struct SolutionState : State {
    bool is_solved;
    char failure_type;
    double GDOP;

    SolutionState();
    SolutionState(int time, 
                  const std::vector<double>& position,
                  const std::vector<double>& velocity,
                  bool is_solved, char failure_type, double GDOP);
};

struct RawMeasurement {
    bool is_present;
    int time;
    unsigned prn_id;
    double L1_range, L2_range;
    double L1_phase, L2_phase;
    unsigned short L1_SNR, L2_SNR;
    unsigned qualflg;

    RawMeasurement();
    RawMeasurement(bool is_present, 
                   int time, unsigned prn_id,
                   double L1_range, double L2_range,
                   double L1_phase, double L2_phase,
                   unsigned short L1_SNR, unsigned short L2_SNR,
                   unsigned qualflg);
};

struct RefinedMeasurement {
    bool is_present;
    int time;
    unsigned prn_id;
    double pseudorange;
    double carrier_phase;
    std::vector<double> gps_position;
    std::vector<double> gps_velocity;

    RefinedMeasurement();
    RefinedMeasurement(bool is_present,
                       int time, unsigned prn_id, double pseudorange, double carrier_phase,
                       const std::vector<double>& gps_position, const std::vector<double>& gps_velocity);
};

struct RawMeasurementGroupped {
    int time;
    std::vector<RawMeasurement> raw_measurements;

    RawMeasurementGroupped();
    RawMeasurementGroupped(int time,
                           const std::vector<RawMeasurement>& raw_measurements);
};

struct RefinedMeasurementGroupped {
    int time;
    std::vector<RefinedMeasurement> refined_measurements;

    RefinedMeasurementGroupped();
    RefinedMeasurementGroupped(int time,
                               const std::vector<RefinedMeasurement>& refined_measurements);
};
