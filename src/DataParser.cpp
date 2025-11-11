#include "DataParser.hpp"

std::vector<State> DataParser::load_grace_fo_gnv_data(const std::string& filename) {
    std::vector<State> true_states;

    std::ifstream file;
    file.open("../data/gnv/" + filename, std::ios::in);

    std::string line;

    unsigned skiprows = 148;
    for (unsigned i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {
        std::stringstream stream(line);

        std::getline(stream, line, ' ');
        unsigned user_time = std::stoi(line);

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

        std::vector<double> position = {xpos, ypos, zpos};
        std::vector<double> velocity = {xvel, yvel, zvel};
        true_states.push_back({user_time, position, velocity});
    }

    file.close();

    return true_states;
}

std::vector<RawMeasurementGroupped> DataParser::load_grace_fo_gps_data(const std::string& filename) {
    std::vector<RawMeasurementGroupped> rmgs;
    std::vector<RawMeasurement> rms(32);
    unsigned user_time_current = 0;
    
    std::ifstream file;
    file.open("../data/gps/" + filename, std::ios::in);

    std::string line;

    unsigned skiprows = 196;
    for (unsigned i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {
        std::stringstream stream(line);

        std::getline(stream, line, ' ');

        unsigned user_time = std::stoi(line);
    
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

        unsigned qualflg = std::stoi(line);

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

        RawMeasurement rm = {user_time, prn_id,
                             L1_range, L2_range,
                             L1_phase, L2_phase,
                             L1_SNR, L2_SNR,
                             qualflg};

        if (user_time_current == 0) user_time_current = user_time;

        if (user_time == user_time_current) {
            unsigned prn_index = prn_id - 1;
            rms[prn_index] = rm;
        } else {
            rmgs.push_back({user_time_current, rms});
            rms.clear();
            user_time_current = user_time;
        }
    }

    file.close();

    return rmgs;
}

void DataParser::load_brdc_data(const std::string& filename, 
                                std::map<unsigned, std::vector<unsigned>>& ts_ephs, 
                                std::map<std::pair<unsigned, unsigned>, Ephemeris>& ephs) {
    std::ifstream file;
    file.open("../data/brdc/" + filename, std::ios::in);

    std::string line;

    unsigned skiprows = 8;
    for (unsigned i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {
    
        unsigned prn_id = std::stoi(line.substr(0, 2));
        double a_f0 = std::stod(line.substr(22, 41));
        double a_f1 = std::stod(line.substr(41, 60));
        double a_f2 = std::stod(line.substr(60, 79));

        std::getline(file, line);

        double C_rs = std::stod(line.substr(22, 41));
        double delta_n = std::stod(line.substr(41, 60));
        double M_0 = std::stod(line.substr(60, 79));

        std::getline(file, line);

        double C_uc = std::stod(line.substr(3, 22));
        double e = std::stod(line.substr(22, 41));
        double C_us = std::stod(line.substr(41, 60));
        double A_sqrt = std::stod(line.substr(60, 79));

        std::getline(file, line);

        unsigned t_oe = std::stod(line.substr(3, 22));
        double C_ic = std::stod(line.substr(22, 41));
        double Omega_0 = std::stod(line.substr(41, 60));
        double C_is = std::stod(line.substr(60, 79));
    
        std::getline(file, line);

        double i_0 = std::stod(line.substr(3, 22));
        double C_rc = std::stod(line.substr(22, 41));
        double omega = std::stod(line.substr(41, 60));
        double Omega_dot = std::stod(line.substr(60, 79));

        std::getline(file, line);

        double IDOT = std::stod(line.substr(3, 22));

        std::getline(file, line);
        std::getline(file, line);

        ts_ephs[prn_id].push_back(t_oe);

        std::pair key = {prn_id, t_oe};
        ephs[key].a_f0 = a_f0;
        ephs[key].a_f0 = a_f0;
        ephs[key].a_f1 = a_f1;
        ephs[key].a_f2 = a_f2;
        ephs[key].M_0 = M_0;
        ephs[key].delta_n = delta_n;
        ephs[key].e = e;
        ephs[key].A_sqrt = A_sqrt;
        ephs[key].Omega_0 = Omega_0;
        ephs[key].i_0 = i_0;
        ephs[key].omega = omega;
        ephs[key].Omega_dot = Omega_dot;
        ephs[key].IDOT = IDOT;
        ephs[key].C_uc = C_uc;
        ephs[key].C_us = C_us;
        ephs[key].C_rc = C_rc;
        ephs[key].C_rs = C_rs;
        ephs[key].C_ic = C_ic;
        ephs[key].C_is = C_is;
    }

    file.close();
}
