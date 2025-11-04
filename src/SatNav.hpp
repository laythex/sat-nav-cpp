#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <limits>

#include "LinAlg.hpp"
#include "GPSHandler.hpp"

class SatNav {

public:
    SatNav(const std::string& gnv_filename, const std::string& gps_filename, 
           const GPSHandler& handler);
    SatNav(const std::string& gnv_filename1, const std::string& gnv_filename2, 
           const std::string& gps_filename1, const std::string& gps_filename2, 
           const GPSHandler& handler);

    void solve(unsigned ti = 0, unsigned tf = 0, bool is_relativistic = true);
    void solve_rel(unsigned ti = 0, unsigned tf = 0);

    void out_errors_norm();
    void out_errors_prs();
    void out_errors_rel();
    void out_number_of_sats();
    void out_snr_over_sol();

private:
    void load_gnv_data(const std::string& gnv_filename, std::map<unsigned, std::vector<double>>& p, std::map<unsigned, std::vector<double>>& v);
    void load_gps_data(const std::string& gps_filename, std::map<std::pair<unsigned, unsigned>, std::vector<double>>& r, std::set<unsigned>& t);

    std::vector<double> correct_raw(unsigned prn_id, unsigned gps_time, const std::vector<double>& measurments, bool is_relativistic = true);
    std::vector<double> calc_pos_from_raw(const std::vector<double>& pseudoranges, const std::vector<std::vector<double>>& gps_positions) const;
    std::vector<double> calc_rel_pos_from_raw(const std::vector<double>& pseudoranges_diff, const std::vector<std::vector<double>>& gps_positions, 
                                              const std::vector<double>& passive_position) const;
    std::vector<bool> find_low_satellites(const std::vector<double>& sol, const std::vector<std::vector<double>>& gps_positions);
    bool is_fading(unsigned prn_id, const std::set<unsigned>::iterator& it_ts);

    double c = 2.99792458e8;
    double earth_rotation_rate = 7.2921151467e-5;
    double C1 = 2.54572778016;
    double C2 = -1.54572778016;

    double GDOP0 = 5;
    double mask_angle = 5;
    unsigned fadeout_time = 20;
    double SNR_threshold = 0;

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
    std::vector<std::pair<unsigned, unsigned>> number_of_sats;
};
