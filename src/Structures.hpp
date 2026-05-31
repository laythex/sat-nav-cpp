#pragma once

#include <vector>

enum class GRACE : char { A = 'A', B = 'B', C = 'C', D = 'D' };

enum class SWARM : char { A = 'A', B = 'B', C = 'C' };

enum class MeasurementStatus {
    DEFAULT,
    OK,
    NOT_PRESENT,
    QUALITY_FLAG,
    CN0,
    FADING,
    COUNT
};

enum class SolutionStatus {
    DEFAULT,
    OK,
    SINGULAR_B1,
    HIGH_GDOP,
    COUNT
};

enum class IntendedError {
    NONE,
    QUALITY_FLAG,
    CN0,
    NEW,
    LOW,
    GPS_CLOCK,
    IONOSPHERIC,
    RELATIVISTIC,
    SAGNAC,
    HATCH_FILTER
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
    std::vector<double> position = std::vector<double>(3, 0.0);
    std::vector<double> velocity = std::vector<double>(3, 0.0);
};

struct SolutionState : State {
    SolutionStatus status = SolutionStatus::DEFAULT;
    double GDOP = 0.0;
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

struct RefinedMeasurement {
    MeasurementStatus status = MeasurementStatus::DEFAULT;
    unsigned prn_id = 0;
    unsigned toa_ASN = 0;
    double tot = 0;
    double pseudorange = 0.0;
    double carrier_phase = 0.0;
    std::vector<double> gps_position = std::vector<double>(3, 0.0);
    std::vector<double> gps_velocity = std::vector<double>(3, 0.0);
};

struct RawMeasurementGroupped {
    unsigned time = 0;
    std::vector<RawMeasurement> raw_measurements = std::vector<RawMeasurement>(32);
};

struct RefinedMeasurementGroupped {
    unsigned time = 0;
    std::vector<RefinedMeasurement> refined_measurements = std::vector<RefinedMeasurement>(32);
};

struct DynamicFilterState {
    unsigned time = 0;
    std::vector<double> x = std::vector<double>(3, 0.0);
    std::vector<double> dx = std::vector<double>(3, 0.0);
    std::vector<double> pas_pos = std::vector<double>(3, 0.0);
    std::vector<double> d_pas_pos = std::vector<double>(3, 0.0);
    std::vector<bool> mask = std::vector<bool>(32, false);
    std::vector<double> xi_m_1 = std::vector<double>(3, 0.0);
    std::vector<double> xi_m_2 = std::vector<double>(3, 0.0);
    std::vector<std::vector<double>> C_rows = std::vector<std::vector<double>>(32);
    std::vector<double> x_est = std::vector<double>(3, 0.0);
};

struct AccelerationMeasurement {
    unsigned time = 0;
    std::vector<double> linear_acceleration = std::vector<double>(3, 0.0);
    std::vector<double> angular_acceleration = std::vector<double>(3, 0.0);
};

constexpr char to_char(GRACE sat) {
    return static_cast<char>(sat);
}

constexpr char to_char(SWARM sat) {
    return static_cast<char>(sat);
}

constexpr std::string to_string(MeasurementStatus status) {
    constexpr const char* messages[] = {
        "Status not assigned",
        "OK",
        "Measurement is not present",
        "Quality flag is bad",
        "CN0 is bad",
        "Satellite is fading"
    };
    return messages[static_cast<int>(status)];
}

constexpr std::string to_string(SolutionStatus status) {
    constexpr const char* messages[] = {
        "Status not assigned",
        "OK",
        "Matrix B1 is singular",
        "GDOP is too high"
    };
    return messages[static_cast<int>(status)];
}

constexpr std::string to_filename(IntendedError status) {
    constexpr const char* messages[] = {
        "none",
        "qualflg",
        "cn0",
        "new",
        "low",
        "gps-clock",
        "ionophere",
        "relativity",
        "sagnac",
        "hatch",
    };
    return messages[static_cast<int>(status)];
}

constexpr std::string to_title(IntendedError status) {
    constexpr const char* messages[] = {
        "none",
        "признак качества",
        "отношение несущей к плотности шума",
        "новые НКА",
        "низкие НКА",
        "ошибка часов НКА",
        "ионосферный вклад",
        "релятивистская поправка",
        "эффект Саньяка",
        "хатч-фильтр",
    };
    return messages[static_cast<int>(status)];
}
