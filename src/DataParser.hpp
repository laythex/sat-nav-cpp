#pragma once

#include <iostream>

#include <fstream>
#include <iomanip>

#include <string>
#include <vector>
#include <map>
#include <set>

#include "Structures.hpp"

class DataParser {
public:
    static std::vector<State> load_grace_fo_gnv_data(const std::string& date, const GRACE_SATS& sat);
    static std::vector<RawMeasurementGroupped> load_grace_fo_gps_data(const std::string& date, const GRACE_SATS& sat);

    static std::vector<State> load_grace_gnv_data(const std::string& date, const GRACE_SATS& sat);
    static std::vector<RawMeasurementGroupped> load_grace_gps_data(const std::string& date, const GRACE_SATS& sat);
    static std::vector<AccelerationMeasurement> load_grace_acc_data(const std::string& date, const GRACE_SATS& sat);
    
    static std::vector<State> load_swarm_nav_data(const std::string& date, const SWARM_SATS& sat);
    static std::vector<RawMeasurementGroupped> load_swarm_gps_data(const std::string& date, const SWARM_SATS& sat);

    static std::vector<std::vector<Ephemeris>> load_brdc_data(const std::string& date);

private:
    static unsigned grace_to_gps_time(unsigned grace_time);
    static unsigned date_to_gps_time(unsigned year, unsigned month, unsigned day, unsigned hours, unsigned minutes, unsigned seconds);

    template <typename T>
    static T bytes_to_T_endian(const char* bytes, T);

    static inline std::string grace_gnv_dir = "../data/grace/gnv/";
    static inline std::string grace_gps_dir = "../data/grace/gps/";
    static inline std::string grace_acc_dir = "../data/grace/acc/";
    static inline std::string swarm_nav_dir = "../data/swarm/nav/";
    static inline std::string swarm_gps_dir = "../data/swarm/gps/";
    static inline std::string gps_brdc_dir = "../data/gps/brdc/";

    static inline unsigned month_day_count[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    static inline unsigned day_seconds_count = 86400;
};

template <typename T>
T DataParser::bytes_to_T_endian(const char* bytes, T) {
    size_t n = sizeof(T);
    T result{};

    for (size_t i = 0; i < n; ++i) {
        reinterpret_cast<char*>(&result)[i] = bytes[n - 1 - i];
    }

    return result;
}
