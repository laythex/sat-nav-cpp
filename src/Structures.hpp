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
    unsigned time;
    std::vector<double> position;
    std::vector<double> velocity;

    State();
    State(unsigned time, 
          const std::vector<double>& position,
          const std::vector<double>& velocity);
};

struct GPSState : State {
    double relativistic_error;

    GPSState();
    GPSState(unsigned time, 
             const std::vector<double>& position,
             const std::vector<double>& velocity,
             double relativistic_error);
};

struct SolutionState : State {
    bool is_solved;
    char failure_type;
    double GDOP;

    SolutionState();
    SolutionState(unsigned time, 
                  const std::vector<double>& position,
                  const std::vector<double>& velocity,
                  bool is_solved, char failure_type, double GDOP);
};

struct RawMeasurement {
    bool is_present;
    unsigned time;
    unsigned prn_id;
    double L1_range, L2_range;
    double L1_phase, L2_phase;
    double L1_SNR, L2_SNR;
    unsigned qualflg;

    RawMeasurement();
    RawMeasurement(bool is_present, 
                   unsigned time, unsigned prn_id,
                   double L1_range, double L2_range,
                   double L1_phase, double L2_phase,
                   double L1_SNR, double L2_SNR,
                   unsigned qualflg);
};

struct RefinedMeasurement {
    bool is_present;
    unsigned time;
    unsigned prn_id;
    double pseudorange;
    double carrier_phase;
    std::vector<double> gps_position;
    std::vector<double> gps_velocity;

    RefinedMeasurement();
    RefinedMeasurement(bool is_present,
                       unsigned time, unsigned prn_id, double pseudorange, double carrier_phase,
                       const std::vector<double>& gps_position, const std::vector<double>& gps_velocity);
};

struct RawMeasurementGroupped {
    unsigned time;
    std::vector<RawMeasurement> raw_measurements;

    RawMeasurementGroupped();
    RawMeasurementGroupped(unsigned time,
                           const std::vector<RawMeasurement>& raw_measurements);
};

struct RefinedMeasurementGroupped {
    unsigned time;
    std::vector<RefinedMeasurement> refined_measurements;

    RefinedMeasurementGroupped();
    RefinedMeasurementGroupped(unsigned time,
                               const std::vector<RefinedMeasurement>& refined_measurements);
};

// Используется в алгоритме динамической фильтрации задачи поиска вектора относительного состояния
struct DynamicFilterState {
    unsigned time;
    std::vector<double> x;
    std::vector<double> dx;
    std::vector<double> pas_pos;
    std::vector<double> d_pas_pos;
    std::vector<bool> mask;
    std::vector<double> xi_m_1;
    std::vector<double> xi_m_2;
    std::vector<std::vector<double>> C_rows;
    std::vector<double> x_est;

    DynamicFilterState();
    DynamicFilterState(unsigned time, 
                       const std::vector<double>& x, const std::vector<double>& dx,
                       const std::vector<double>& pas_pos, const std::vector<double>& d_pas_pos,
                       const std::vector<bool>& mask, const std::vector<double>& xi_m_1, const std::vector<double>& xi_m_2,
                       const std::vector<std::vector<double>>& C_rows, const std::vector<double>& x_est);
};

struct AccelerationMeasurement {
    unsigned time;
    std::vector<double> linear_acceleration;
    std::vector<double> angular_acceleration;

    AccelerationMeasurement();
    AccelerationMeasurement(unsigned time,
                            const std::vector<double>& linear_acceleration,
                            const std::vector<double>& angular_acceleration);
};

enum class GRACE_SATS {
    A, B, C, D
};

enum class SWARM_SATS {
    A, B, C
};

constexpr char satToChar(GRACE_SATS sat) {
    switch (sat) {
        case GRACE_SATS::A:   return 'A';
        case GRACE_SATS::B:   return 'B';
        case GRACE_SATS::C:   return 'C';
        case GRACE_SATS::D:   return 'D';
        default: return '0';
    }
}

constexpr char satToChar(SWARM_SATS sat) {
    switch (sat) {
        case SWARM_SATS::A:   return 'A';
        case SWARM_SATS::B:   return 'B';
        case SWARM_SATS::C:   return 'C';
        default: return '0';
    }
}

