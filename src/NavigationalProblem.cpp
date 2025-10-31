#include "NavigationalProblem.hpp"

NavigationalProblem::NavigationalProblem(const std::string& gnv_filename, const std::string& gps_filename, const GPSHandler& handler) : handler(handler) {
    load_gnv_data(gnv_filename);
    load_gps_data(gps_filename);
}

void NavigationalProblem::load_gnv_data(const std::string& gnv_filename) {
    std::ifstream gnv_file;
    gnv_file.open(gnv_filename, std::ios::in);

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

        pos[gps_time] = {xpos, ypos, zpos};
        vel[gps_time] = {xvel, yvel, zvel};
    }

    gnv_file.close();
}

void NavigationalProblem::load_gps_data(const std::string& gps_filename) {
    std::ifstream gps_file;
    gps_file.open(gps_filename, std::ios::in);

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
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        double L1_range = std::stod(line);
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        double L2_range = std::stod(line);

        ts.insert(gps_time);
        raw[{prn_id, gps_time}] = {L1_range, L2_range};
    }

    gps_file.close();
}

void NavigationalProblem::solve(unsigned ti, unsigned tf) {

    std::set<unsigned> ts_range;
    auto it_i = (ti > 0) ? ts.lower_bound(ti) : ts.begin();
    auto it_f = (tf > 0) ? ts.upper_bound(tf) : ts.end();
    for (auto it_ts = it_i; it_ts != it_f; it_ts++) {
        ts_range.insert(*it_ts);
    }
    ts = ts_range;

    for (auto it_ts = ts.begin(); it_ts != ts.end(); it_ts++) {
        unsigned t = *it_ts;

        std::vector<double> pseudoranges;
        std::vector<std::vector<double>> gps_positions;

        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            auto it_raw = raw.find({prn_id, t});
            if (it_raw == raw.end()) continue;    
            
            std::vector<double> measurments = it_raw->second;
            std::vector<double> corrected = correct(prn_id, t, measurments);
            
            double pseudorange = corrected[3];
            std::vector<double> gps_pos(corrected.begin(), corrected.begin() + 3);

            pseudoranges.push_back(pseudorange);
            gps_positions.push_back(gps_pos);
            
            err_prs[{prn_id, t}] = std::abs(pseudorange - abs(gps_pos - pos[t]));
        }

        try {
            err_sol[t] = iterative(pseudoranges, gps_positions) - pos[t];
        } catch (const std::runtime_error&) {
            ts.erase(t);
            it_ts--;
        }

    }
}

std::vector<double> NavigationalProblem::correct(unsigned prn_id, unsigned gps_time, const std::vector<double>& measurments) {
    double L1_range = measurments[0];
    double L2_range = measurments[1];

    double iono_free = C1 * L1_range + C2 * L2_range;
    
    /*  Эффект                     t       t * v    t * c   
    (1) Задержка распространения - 80 мс - 240 км -  ---   
    (2) Ошибка часов НКА         - 3 мс  - 10 м   - 1000 км 
    (3) Релятивизм               - 1 мкс - 3 мм   - 300 м  

    При расчете ошибки часов учитываем только (1)
    При расчете эфемерид учитываем (1), (2) 
    При корректировке псевдодальности учитываем (1), (2), (3) */

    double delay = 0;

    double propagation_delay = iono_free / c;
    delay += propagation_delay;

    double clock_error = handler.get_clock_error(prn_id, gps_time - delay);
    delay += clock_error;

    std::vector<double> state = handler.get_state(prn_id, gps_time - delay);

    double relativity_error = state[6];
    delay += relativity_error;

    double pseudorange = delay * c;
    std::vector<double> gps_pos(state.begin(), state.begin() + 3);

    double phi = delay * earth_rotation_rate;
    Matrix rot = rotation(-phi, 'z');
    gps_pos = rot * gps_pos;

    std::vector<double> corrected = {gps_pos[0], gps_pos[1], gps_pos[2], pseudorange};
    return corrected;
}

std::vector<double> NavigationalProblem::iterative(const std::vector<double>& pseudoranges, const std::vector<std::vector<double>>& gps_positions) const {

    std::vector<double> x0 = {0.0, 0.0, 0.0};
    double c_tau = 0.0;
    double delta, eps = 1.0;
    
    unsigned n = pseudoranges.size();

    std::vector<double> U(n);
    Matrix B(n, 4, 1.0);
    std::vector<double> dX(4);
    std::vector<double> dX_X(3);
    double dX_c_tau;

    std::vector<double> DX(n);
    double D;

    while (true) {
        
        for (unsigned i = 0; i < n; i++) {
            DX = gps_positions[i] - x0;
            D = abs(DX);

            U[i] = pseudoranges[i] - D;
            for (unsigned j = 0; j < 3; j++) {
                B.at(i, j) = DX[j] / D;
            }
        }

        dX = (B.transpose() * B).inverse() * B.transpose() * U * (-1.0);

        dX_X = std::vector<double>(dX.begin(), dX.begin() + 3);
        dX_c_tau = dX[3];
    
        delta = abs(dX_X) + std::abs(dX_c_tau - c_tau);
        if (delta < eps) {
            break;
        }

        x0 = x0 + dX_X;
        c_tau = dX_c_tau;
    }

    return x0;
}

void NavigationalProblem::errors_norm(const std::string& err_filename) {
    std::ofstream err_file;
    err_file.open(err_filename, std::fstream::out);

    err_file << "Модуль ошибки" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;

    for (const unsigned& t : ts) {
        double error = abs(err_sol[t]);
        err_file << t << '\t' << error << std::endl;
    }

    err_file.close();
}

void NavigationalProblem::errors_prs(const std::string& err_filename) {
    std::ofstream err_file;
    err_file.open(err_filename, std::fstream::out);

    err_file << "Модуль ошибки псевдодальности" << '\t' << "Время, с" << '\t' << "Ошибка, м" << std::endl;

    for (const unsigned& t : ts) {
        err_file << t << '\t';

        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            auto it = err_prs.find({prn_id, t});

            if (it != err_prs.end()) {
                err_file << err_prs[{prn_id, t}];
            }

            err_file << '\t';
        }

        err_file << '\n';
    }

    err_file.close();
}