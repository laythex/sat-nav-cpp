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
    
    std::vector<std::vector<Ephemeris>> load_brdc_data(const std::string& filename);
};