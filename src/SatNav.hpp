#pragma once

#include <iostream>

#include <vector>
#include <map>
#include <set>

#include <stdexcept>
#include <algorithm>

#include "LinAlg.hpp"
#include "GPSHandler.hpp"
#include "DataParser.hpp"

class SatNav {

public:
    SatNav(const std::string& gnv_filename, const std::string& gps_filename, 
           const GPSHandler& handler);
    SatNav(const std::string& gnv_filename1, const std::string& gnv_filename2, 
           const std::string& gps_filename1, const std::string& gps_filename2, 
           const GPSHandler& handler);

    void solve(unsigned ti = 0, unsigned tf = 0);
    void solve_rel(unsigned ti = 0, unsigned tf = 0);

    void out_error_norm(); // Все out в отдельный файл?
    void out_error_prs();
    void out_number_of_sats();
    void out_is_solved();
    void out_error_by_type(unsigned et);

private:
    std::vector<double> correct_raw(unsigned prn_id, unsigned gps_time, const std::vector<double>& L_ranges);
    std::vector<double> calc_pos_from_raw(const std::vector<double>& pseudoranges, const std::vector<std::vector<double>>& gps_positions) const;
    std::vector<double> calc_rel_pos_from_raw(const std::vector<double>& pseudoranges_diff, const std::vector<std::vector<double>>& gps_positions, 
                                              const std::vector<double>& passive_position) const;
    std::vector<bool> mask_low_satellites(const std::vector<double>& sol, const std::vector<std::vector<double>>& gps_positions);
    bool is_fading(unsigned prn_id, const std::set<unsigned>::iterator& it_ts);

    double c = 2.99792458e8;
    double earth_rotation_rate = 7.2921151467e-5;
    double gamma = 1.64694444444;
    double C1 = 2.54572778016;
    double C2 = -1.54572778016;

    double GDOP0 = 5;
    double mask_angle = 5;
    unsigned fadeout_time = 20;
    double SNR_threshold = 30;
    double hatch_constant = 0.1;

    GPSHandler handler;

    std::map<unsigned, std::vector<double>> pos;
    std::map<unsigned, std::vector<double>> vel;
    std::map<unsigned, std::vector<double>> pos2;
    std::map<unsigned, std::vector<double>> vel2;

    std::set<unsigned> ts_raw;
    std::map<std::pair<unsigned, unsigned>, std::vector<double>> raw;
    std::map<std::pair<unsigned, unsigned>, std::vector<double>> raw2;

    std::set<unsigned> ts_sol;
    std::map<std::pair<unsigned, unsigned>, double> err_prs;
    std::map<unsigned, std::vector<double>> err_sol;
    std::map<unsigned, unsigned> number_of_sats;
    std::map<unsigned, bool> is_solved;
    
    std::vector<std::string> error_names = {"iono", "rel", "fade", "low", "snr", "hatch"};
    std::vector<std::string> error_descr = {"Ионосферная комбинация",
                                            "Релятивистская поправка",
                                            "Появляющиеся/исчезающие спутники, (t > " + std::to_string(fadeout_time) + " с)",
                                            "Низкие спутники, (высота > " + std::to_string(mask_angle) + " град)",
                                            "Сигнал/шум (SNR > " + std::to_string(SNR_threshold) + ")",
                                            "Hatch-фильтр"};
    unsigned error_type = -1;
};
