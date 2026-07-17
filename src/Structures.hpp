#pragma once

#include <Eigen/Dense>
#include <format>
#include <vector>

#include "Constants.hpp"

enum class GRACE : char { A = 'A', B = 'B', C = 'C', D = 'D' };

enum class SWARM : char { A = 'A', B = 'B', C = 'C' };

enum class MeasurementStatus { DEFAULT, OK, NOT_PRESENT, QUALITY_FLAG, CN0, FADE_IN };

enum class SolutionStatus { DEFAULT, OK, SINGULAR_B1, HIGH_GDOP };

enum class SolutionStatusRel { DEFAULT, OK, NO_STANDALONE, INITIALIZING, MODEL_STEP, SINGULAR_C1 };

enum class IntendedError {
    NONE,
    QUALITY_FLAG,
    NOISY,
    GPS_CLOCK,
    IONOSPHERIC,
    RELATIVISTIC,
    SAGNAC,
    HATCH,
    GDOP,
    DELTA3,
    GAMMA,
    KALMAN
};

struct Date {
    unsigned year;
    unsigned month;
    unsigned day;
};

struct Ephemeris {
    unsigned prn_id;
    int t_oe;

    double a_f0, a_f1, a_f2;
    double M_0, delta_n, e, A_sqrt, Omega_0, i_0, omega, Omega_dot, IDOT;
    double C_uc, C_us, C_rc, C_rs, C_ic, C_is;
};

struct State {
    unsigned time = 0;
    Eigen::Vector3d position = Eigen::Vector3d::Zero();
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
};

struct DiscreteState {
    unsigned time = 0;
    Eigen::Vector3d position = Eigen::Vector3d::Zero();
    Eigen::Vector3d increment = Eigen::Vector3d::Zero();
};

struct RawMeasurement {
    bool is_present = false;
    unsigned time = 0;
    unsigned prn_id = 0;
    double L1_range = 0.0;
    double L2_range = 0.0;
    double L1_phase = 0.0;
    double L2_phase = 0.0;
    double L1_CN0 = 0.0;
    double L2_CN0 = 0.0;
    unsigned qualflg = 0;
};

struct RawMeasurementGroupped {
    unsigned time = 0;
    std::vector<RawMeasurement> measurements = std::vector<RawMeasurement>(32);
};

struct RefinedMeasurement {
    MeasurementStatus status = MeasurementStatus::DEFAULT;
    unsigned prn_id = 0;
    unsigned toa_ASN = 0;
    double tot = 0.0;
    double pr_refined = 0.0;
    double cp_refined = 0.0;
    State gps_state;
};

struct RefinedMeasurementGroupped {
    unsigned time = 0;
    std::vector<RefinedMeasurement> measurements = std::vector<RefinedMeasurement>(32);
};

struct AccelerationMeasurement {
    unsigned time = 0;
    std::vector<double> linear_acceleration = std::vector<double>(3, 0.0);
    std::vector<double> angular_acceleration = std::vector<double>(3, 0.0);
};

struct SolutionState : State {
    SolutionStatus status = SolutionStatus::DEFAULT;
    double GDOP = 0.0;
    double dt_ASN_c = 0.0;
    double residual = 0.0;

    std::vector<bool> gps_mask = std::vector<bool>(32, false);
    std::vector<State> gps_states = std::vector<State>(32);
    std::vector<double> prs_L1 = std::vector<double>(32, 0.0);
    std::vector<double> cps_L1 = std::vector<double>(32, 0.0);
};

struct SolutionStateRel : State {
    SolutionStatusRel status = SolutionStatusRel::DEFAULT;
};

struct StandaloneData {
    unsigned time = 0;
    SolutionState pas;
    SolutionState act;

    bool is_fully_solved() const;
    bool is_fully_present_at(size_t prn_id) const;
};

constexpr char to_char(GRACE sat) { return static_cast<char>(sat); }

constexpr char to_char(SWARM sat) { return static_cast<char>(sat); }

constexpr std::string to_string(MeasurementStatus status) {
    constexpr const char* messages[] = {"Status not assigned", "OK",      "Measurement not present",
                                        "Bad quality flag",    "Bad CN0", "Fading satellite"};
    return messages[static_cast<int>(status)];
}

constexpr std::string to_string(SolutionStatus status) {
    constexpr const char* messages[] = {"Status not assigned", "OK", "Singular B1", "High GDOP"};
    return messages[static_cast<int>(status)];
}

constexpr std::string to_string(SolutionStatusRel status) {
    constexpr const char* messages[] = {"Status not assigned",    "OK",         "Standalone not solved",
                                        "Filter is initializing", "Model step", "Singular C1"};
    return messages[static_cast<int>(status)];
}

constexpr std::string to_filename(IntendedError status) {
    constexpr const char* messages[] = {"none",   "qualflg", "noisy", "gps-clock", "ionospheric", "relativistic",
                                        "sagnac", "hatch",   "gdop",  "delta3",    "gamma",       "kalman"};
    return messages[static_cast<int>(status)];
}

constexpr std::string to_title(IntendedError status) {
    constexpr const char* messages[] = {"none",
                                        "признак качества",
                                        "шумные измерения",
                                        "ошибка часов НКА",
                                        "ионосферный вклад",
                                        "релятивистская поправка",
                                        "эффект Саньяка",
                                        "хатч-фильтр",
                                        "GDOP",
                                        "delta3",
                                        "gamma",
                                        "kalman"};
    return messages[static_cast<int>(status)];
}

constexpr std::string to_string(const Date& date) {
    return std::format("{:04}-{:02}-{:02}", date.year, date.month, date.day);
}
