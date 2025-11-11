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
    std::vector<State> load_grace_fo_gnv_data(const std::string& gnv_filename);

    std::vector<RawMeasurementGroupped> load_grace_fo_gps_data(const std::string& gps_filename);
    
    void load_brdc_data(const std::string& filename, // Сделать констом
                        std::map<unsigned, std::vector<unsigned>>& ts_ephs, 
                        std::map<std::pair<unsigned, unsigned>, Ephemeris>& ephs);
};