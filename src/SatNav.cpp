#include "SatNav.hpp"

SatNav::SatNav(const std::string& gnv_filename, const std::string& gps_filename, 
               const GPSHandler& handler) : handler(handler) {
    
    load_gnv_data(gnv_filename, pos, vel);
    load_gps_data(gps_filename, raw, ts_raw);
}

SatNav::SatNav(const std::string& gnv_filename1, const std::string& gnv_filename2, 
               const std::string& gps_filename1, const std::string& gps_filename2, 
               const GPSHandler& handler) : handler(handler) {
    
    load_gnv_data(gnv_filename1, pos, vel);
    load_gnv_data(gnv_filename2, pos2, vel2);

    std::set<unsigned> t1, t2;

    load_gps_data(gps_filename1, raw, t1);
    load_gps_data(gps_filename2, raw2, t2);

    for (unsigned t : t1) {
        if (t2.find(t) != t2.end()) ts_raw.insert(t);
    }
}

void SatNav::load_gnv_data(const std::string& gnv_filename, std::map<unsigned, std::vector<double>>& p, std::map<unsigned, std::vector<double>>& v) {
    std::ifstream gnv_file;
    gnv_file.open("../data/gnv/" + gnv_filename, std::ios::in);

    std::string line;

    unsigned skiprows = 148;
    for (unsigned i = 0; i < skiprows; ++i) {
        std::getline(gnv_file, line);
    }

    while (std::getline(gnv_file, line)) {
        std::stringstream stream(line);

        std::getline(stream, line, ' ');
        unsigned gps_time = std::stoi(line);

        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        std::getline(stream, line, ' ');
        double xpos = std::stod(line);
        std::getline(stream, line, ' ');
        double ypos = std::stod(line);
        std::getline(stream, line, ' ');
        double zpos = std::stod(line);
    
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        
        std::getline(stream, line, ' ');
        double xvel = std::stod(line);
        std::getline(stream, line, ' ');
        double yvel = std::stod(line);
        std::getline(stream, line, ' ');
        double zvel = std::stod(line);

        p[gps_time] = {xpos, ypos, zpos};
        v[gps_time] = {xvel, yvel, zvel};
    }

    gnv_file.close();
}

void SatNav::load_gps_data(const std::string& gps_filename, std::map<std::pair<unsigned, unsigned>, std::vector<double>>& r, std::set<unsigned>& t) {
    std::ifstream gps_file;
    gps_file.open("../data/gps/" + gps_filename, std::ios::in);

    std::string line;

    unsigned skiprows = 196;
    for (unsigned i = 0; i < skiprows; ++i) {
        std::getline(gps_file, line);
    }

    while (std::getline(gps_file, line)) {
        std::stringstream stream(line);

        std::getline(stream, line, ' ');

        unsigned gps_time = std::stoi(line);
    
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        unsigned prn_id = std::stoi(line);

        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        double qualfig = std::stod(line);

        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        double L1_range = std::stod(line);

        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        double L2_range = std::stod(line);

        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        double L1_phase = std::stod(line);

        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        double L2_phase = std::stod(line);

        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        double L1_SNR = std::stod(line);
        
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        double L2_SNR = std::stod(line);

        t.insert(gps_time);
        r[{prn_id, gps_time}] = {L1_range, L2_range, L1_phase, L2_phase, L1_SNR, L2_SNR, qualfig};
    }

    gps_file.close();
}

void SatNav::solve(unsigned ti, unsigned tf) {
    std::vector<double> hatch_L1_range(33);
    std::vector<double> hatch_L2_range(33);
    std::vector<double> hatch_L1_phase(33);
    std::vector<double> hatch_L2_phase(33);
    std::vector<bool> hatch_reset(33, true);

    for (auto it_ts = ts_raw.begin(); it_ts != ts_raw.end(); it_ts++) {
        unsigned t = *it_ts;

        if (ti > 0 || tf > 0) {
            if (t < ti || t >= tf) {
                continue;
            }
        }

        std::vector<double> PR;
        std::vector<std::vector<double>> X;

        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            auto it_raw = raw.find({prn_id, t});
            
            if (it_raw == raw.end()) { hatch_reset[prn_id] = true; continue; }
            if (is_fading(prn_id, it_ts) && error_type != 2) { hatch_reset[prn_id] = true; continue; }

            std::vector<double> measurments = it_raw->second;
            double L1_range = measurments[0];
            double L2_range = measurments[1];
            double L1_phase = measurments[2];
            double L2_phase = measurments[3];
            double L1_SNR = measurments[4];
            double L2_SNR = measurments[5];
            double qualfig = measurments[6];

            if (qualfig != 0) { hatch_reset[prn_id] = true; continue; }
            if ((L1_SNR < SNR_threshold || L2_SNR < SNR_threshold) && error_type != 4) { hatch_reset[prn_id] = true; continue; }

            if (0 && error_type != 5) {
                if (!hatch_reset[prn_id]) {
                    L1_range = hatch_filter(L1_range, hatch_L1_range[prn_id], 
                                            L1_phase, hatch_L1_phase[prn_id]);
                    L2_range = hatch_filter(L2_range, hatch_L2_range[prn_id], 
                                            L2_phase, hatch_L2_phase[prn_id]);
                }
                
                hatch_L1_range[prn_id] = L1_range;
                hatch_L2_range[prn_id] = L2_range;
                hatch_L1_phase[prn_id] = L1_phase;
                hatch_L2_phase[prn_id] = L2_phase;
                hatch_reset[prn_id] = false;
            }

            std::vector<double> corrected = correct_raw(prn_id, t, {L1_range, L2_range});
            double pseudorange = corrected[3];
            std::vector<double> gps_pos(corrected.begin(), corrected.begin() + 3);

            PR.push_back(pseudorange);
            X.push_back(gps_pos);
            
            err_prs[{prn_id, t}] = std::abs(pseudorange - abs(gps_pos - pos[t]));
        }

        number_of_sats.push_back({t, PR.size()}); // использовать ts_raw

        try {
            std::vector<double> sol = calc_pos_from_raw(PR, X);

            if (error_type != 3) {
                std::vector<bool> low_mask = find_low_satellites(sol, X);
                if (std::find(low_mask.begin(), low_mask.end(), true) != low_mask.end()) {
                    PR = mask(PR, low_mask);
                    X = mask(X, low_mask);
                    sol = calc_pos_from_raw(PR, X);
                }
            }

            err_sol[t] = sol - pos[t];
            ts_sol.insert(t);

        } catch (const std::runtime_error& e) {
            continue;
        }
    }
}

void SatNav::solve_rel(unsigned ti, unsigned tf) {

    for (auto it_ts = ts_raw.begin(); it_ts != ts_raw.end(); it_ts++) {
        unsigned t = *it_ts;

        if (ti > 0 || tf > 0) {
            if (t < ti || t >= tf) {
                continue;
            }
        }

        std::vector<double> pseudoranges;
        std::vector<double> pseudoranges_diff;
        std::vector<std::vector<double>> gps_positions;

        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            auto it_raw1 = raw.find({prn_id, t});
            auto it_raw2 = raw2.find({prn_id, t});

            if (it_raw1 == raw.end() || it_raw2 == raw.end()) continue;
            // if (is_fading(prn_id, it_ts)) continue;

            std::vector<double> measurments1 = it_raw1->second;
            std::vector<double> corrected = correct_raw(prn_id, t, measurments1);
            double pseudorange = corrected[3];
            std::vector<double> gps_pos(corrected.begin(), corrected.begin() + 3);

            pseudoranges.push_back(pseudorange);
            gps_positions.push_back(gps_pos);

            std::vector<double> measurments2 = it_raw2->second;

            pseudoranges_diff.push_back(measurments1[0] - measurments2[0]);
        }

        try {
            std::vector<double> sol_passive = calc_pos_from_raw(pseudoranges, gps_positions);

            std::vector<bool> low_mask = find_low_satellites(sol_passive, gps_positions);
            pseudoranges = mask(pseudoranges, low_mask);
            gps_positions = mask(gps_positions, low_mask);
            pseudoranges_diff = mask(pseudoranges_diff, low_mask);

            sol_passive = calc_pos_from_raw(pseudoranges, gps_positions);
            std::vector<double> sol = calc_rel_pos_from_raw(pseudoranges_diff, gps_positions, sol_passive);
            err_sol[t] = sol - pos[t];
            ts_sol.insert(t);

        } catch (const std::runtime_error&) {
            continue;
        }
    }
}

std::vector<double> SatNav::correct_raw(unsigned prn_id, unsigned gps_time, const std::vector<double>& L_ranges) {
    
    /*  Эффект                     t       t * v    t * c   
    (1) Задержка распространения - 80 мс - 240 км -  ---
    (2) Ошибка часов НКА         - 3 мс  - 10 м   - 1000 км
    (3) Релятивизм               - 20 нс - 60 мкм - 6 м

    При расчете ошибки часов учитываем только (1)
    При расчете эфемерид учитываем (1), (2) 
    При корректировке псевдодальности учитываем (1), (2), (3) */

    std::vector<double> gps_pos;
    std::vector<double> corrected_measurments(2);

    for (unsigned i = 0; i < 2; i++) {

        double delay = 0;

        double propagation_delay = L_ranges[i] / c;
        delay += propagation_delay;

        double clock_error = handler.get_clock_error(prn_id, gps_time - delay);
        delay += clock_error;

        std::vector<double> state = handler.get_state(prn_id, gps_time - delay);

        if (error_type != 1) {
            double relativity_error = state[6];
            delay += relativity_error;
        }

        corrected_measurments[i] = delay * c;

        if (i == 0) gps_pos = std::vector<double>(state.begin(), state.begin() + 3);
    }

    double phi = corrected_measurments[0] / c * earth_rotation_rate;
    Matrix rot = rotation(-phi, 'z');
    gps_pos = rot * gps_pos;
    
    double pseudorange = corrected_measurments[0];
    if (error_type != 0) {
        pseudorange = corrected_measurments[0] * C1 + corrected_measurments[1] * C2;
    }

    return {gps_pos[0], gps_pos[1], gps_pos[2], pseudorange};
}

std::vector<double> SatNav::calc_pos_from_raw(const std::vector<double>& PR, const std::vector<std::vector<double>>& X) const {

    std::vector<double> x0 = {0.0, 0.0, 0.0};
    double c_tau = 0.0;
    double eps = 1.0;
    
    unsigned n = PR.size();

    std::vector<double> U(n);
    Matrix B(n, 4, 1.0);

    while (true) {
        for (unsigned i = 0; i < n; i++) {
            std::vector<double> DX = X[i] - x0;
            double D = abs(DX);

            U[i] = PR[i] - D;
            for (unsigned j = 0; j < 3; j++) {
                B.at(i, j) = DX[j] / D;
            }
        }

        Matrix B1 = (B.transpose() * B).inverse();

        double GDOP = sqrt(B1.trace());
        if (GDOP > GDOP0) {
            throw std::runtime_error("GDOP is too high");
        }

        std::vector<double> dX = B1 * B.transpose() * U * (-1.0);

        std::vector<double> dX_X = std::vector<double>(dX.begin(), dX.begin() + 3);
        double dX_c_tau = dX[3];
    
        double delta = abs(dX_X) + std::abs(dX_c_tau - c_tau);
        if (delta < eps) {
            break;
        }

        x0 = x0 + dX_X;
        c_tau = dX_c_tau;
    }

    return x0;
}

std::vector<double> SatNav::calc_rel_pos_from_raw(const std::vector<double>& pseudoranges_diff, const std::vector<std::vector<double>>& gps_positions, 
                                                  const std::vector<double>& passive_position) const {
    unsigned n = pseudoranges_diff.size();

    std::vector<double> U(n);
    Matrix B(n, 3);
        
    for (unsigned i = 0; i < n; i++) {
        unsigned j = (i < n - 1) ? i + 1 : 0;

        double Di = abs(gps_positions[i] - passive_position);
        double Dj = abs(gps_positions[j] - passive_position);

        U[i] = (pseudoranges_diff[i] - pseudoranges_diff[j]);// -( delta[i] - delta[j]);
        
        for (unsigned k = 0; k < 3; k++) {
            B.at(i, k) = (gps_positions[i][k] - passive_position[k]) / Di -
                         (gps_positions[j][k] - passive_position[k]) / Dj;
        }

    }

    std::vector<double> dX = (B.transpose() * B).inverse() * B.transpose() * U * (-1.0);

    return dX;
}

std::vector<bool> SatNav::find_low_satellites(const std::vector<double>& sol, const std::vector<std::vector<double>>& gps_positions) {
    std::vector<bool> mask_indices(gps_positions.size(), true);

    for (unsigned i = 0; i < gps_positions.size(); i++) {
        std::vector<double> gps_rel = gps_positions[i] - sol;
        double zenith_angle = angle_between(sol, gps_rel) * 180.0 * M_1_PI;
        if (90.0 - zenith_angle < mask_angle) {
            mask_indices[i] = false;
            std::cout << "aaa" << std::endl;
        }
    }

    return mask_indices;
}

bool SatNav::is_fading(unsigned prn_id, const std::set<unsigned>::iterator& it_ts) {
    unsigned t0 = *it_ts;

    auto it_ts_fwd = it_ts;
    while (true) {
        it_ts_fwd++;
        if (it_ts_fwd == ts_raw.end()) return true;

        unsigned t = *it_ts_fwd;
        if (t - t0 > fadeout_time) break;

        auto it_raw = raw.find({prn_id, t});
        if (it_raw == raw.end()) return true;
    }

    auto it_ts_bwd = it_ts;
    while (true) {
        if (it_ts_bwd == ts_raw.begin()) return true;
        it_ts_bwd--;

        unsigned t = *it_ts_bwd;
        if (t0 - t > fadeout_time) break;

        auto it_raw = raw.find({prn_id, t});
        if (it_raw == raw.end()) return true;
    }

    return false;
}

double SatNav::hatch_filter(double range, double range_prev, double phase, double phase_prev) {
    return hatch_constant * range + (1 - hatch_constant) * (range_prev + phase - phase_prev);
}

void SatNav::out_error_norm() {
    std::ofstream out_file;
    out_file.open("../results/errors-norm.csv", std::fstream::out);

    out_file << "Модуль ошибки" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;
    out_file << 0 << '\t' << 5 << std::endl;
    
    for (const unsigned& t : ts_sol) {
        out_file << t << ',' << abs(err_sol[t]) << std::endl;
    }

    out_file.close();
}

void SatNav::out_error_prs() {
    std::ofstream out_file;
    out_file.open("../results/errors-prs.csv", std::fstream::out);

    out_file << "Модуль ошибки псевдодальности" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;
    out_file << 0 << '\t' << 5 << std::endl;

    for (const unsigned& t : ts_sol) {
        out_file << t << ',';

        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            auto it = err_prs.find({prn_id, t});

            if (it != err_prs.end()) {
                out_file << err_prs[{prn_id, t}];
            }

            out_file << ',';
        }

        out_file << '\n';
    }

    out_file.close();
}

void SatNav::out_number_of_sats() {
    std::ofstream out_file;
    out_file.open("../results/number-of-sats.csv", std::fstream::out);

    out_file << "Количество видимых спутников" << '\t' << "Время, с" << '\t' << "Кол-во" << std::endl;
    out_file << 0 << '\t' << 15 << std::endl;

    for (auto& num : number_of_sats) {
        out_file << num.first << ',' << num.second << std::endl;
    }

    out_file.close();
}

void SatNav::out_error_by_type(unsigned et) {
    std::ofstream out_file;
    out_file.open("../results/errors-" + error_names[et] + ".csv", std::fstream::out);

    out_file << error_descr[et] << '\t' << "Время, с" << '\t' << "Вклад, м" << std::endl;
    out_file << 0 << '\t' << 0 << std::endl;

    error_type = et;
    solve(0, 0);
    std::map<unsigned, std::vector<double>> err_false = err_sol;

    error_type = -1;
    solve();

    for (const unsigned& t : ts_sol) {
        double error = abs(err_false[t] - err_sol[t]);
        out_file << t << ',' << error << std::endl;
    }

    out_file.close();
}
