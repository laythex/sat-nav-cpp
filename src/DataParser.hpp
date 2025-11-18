#pragma once

#include <iostream>

#include <fstream>
#include <iomanip>

#include <string>
#include <vector>
#include <map>
#include <set>

#include "Structures.hpp"

namespace DataParser {
    std::vector<State> load_grace_fo_gnv_data(const std::string& filename);

    std::vector<RawMeasurementGroupped> load_grace_fo_gps_data(const std::string& filename);

    std::vector<State> load_grace_gnv_data(const std::string& filename);

    std::vector<RawMeasurementGroupped> load_grace_gps_data(const std::string& filename);
    
    void load_brdc_data(const std::string& filename, // Сделать констом
                        std::map<unsigned, std::vector<unsigned>>& ts_ephs, 
                        std::map<std::pair<unsigned, unsigned>, Ephemeris>& ephs);

    // сделать template
    double char_to_double_end(const char* ch);
    double char_to_int_end(const char* ch);
};