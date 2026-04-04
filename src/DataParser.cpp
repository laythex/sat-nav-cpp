#include "DataParser.hpp"

template <typename T>
T bytes_to_T_endian(const char* bytes, T t) {
    size_t n = sizeof(T);
    char res_bytes[n];

    for (size_t i = 0; i < n; i++) {
        res_bytes[n - i - 1] = bytes[i];
    }

    return *reinterpret_cast<T*>(res_bytes);
}

std::vector<State> DataParser::load_grace_fo_gnv_data(const std::string& filename) {
    std::vector<State> true_states;

    std::ifstream file;
    file.open("../data/grace/gnv/" + filename, std::ios::in);

    std::string line;

    unsigned skiprows = 148;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {
        std::stringstream stream(line);

        std::getline(stream, line, ' ');
        int gps_time = std::stoi(line);

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
        true_states.push_back({gps_time, position, velocity});
    }

    file.close();

    return true_states;
}

std::vector<RawMeasurementGroupped> DataParser::load_grace_fo_gps_data(const std::string& filename) {
    std::vector<RawMeasurementGroupped> raw_mgs;
    std::vector<RawMeasurement> raw_ms(32);
    int gps_time_current = 0;
    
    std::ifstream file;
    file.open("../data/grace/gps/" + filename, std::ios::in);

    std::string line;

    unsigned skiprows = 196;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {
        std::stringstream stream(line);

        std::getline(stream, line, ' ');

        int gps_time = std::stoi(line);
    
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

        RawMeasurement raw_m = {true,
                                gps_time, prn_id,
                                L1_range, L2_range,
                                L1_phase, L2_phase,
                                L1_SNR, L2_SNR,
                                qualflg};

        if (gps_time_current == 0) gps_time_current = gps_time;

        if (gps_time != gps_time_current) {
            raw_mgs.push_back({gps_time_current, raw_ms});
            
            raw_ms.clear();
            raw_ms.resize(32);

            gps_time_current = gps_time;
        }

        unsigned prn_index = prn_id - 1;
        raw_ms[prn_index] = raw_m;
    }

    file.close();

    return raw_mgs;
}

#include "LinAlg.hpp"

std::vector<State> DataParser::load_grace_gnv_data(const std::string& filename) {
    std::vector<State> true_states;

    std::ifstream file;
    file.open("../data/grace/gnv/" + filename, std::ios::binary);

    std::string line;
    unsigned skiprows = 26;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    std::ofstream f;
    f.open("../gnv.txt");

    char record[103];
    while (!file.eof()) {
        file.read(record, 4);
        int gps_time = bytes_to_T_endian(record, int());

        if (static_cast<int>(record[0]) == 62) break;

        file.read(record, 2);

        file.read(record, 8);
        double xpos = bytes_to_T_endian(record, double());
        file.read(record, 8);
        double ypos = bytes_to_T_endian(record, double());
        file.read(record, 8);
        double zpos = bytes_to_T_endian(record, double());

        file.read(record, 24);

        file.read(record, 8);
        double xvel = bytes_to_T_endian(record, double());
        file.read(record, 8);
        double yvel = bytes_to_T_endian(record, double());
        file.read(record, 8);
        double zvel = bytes_to_T_endian(record, double());

        file.read(record, 25);

        std::vector<double> position = {xpos, ypos, zpos};
        std::vector<double> velocity = {xvel, yvel, zvel};
        true_states.push_back({gps_time, position, velocity});

        f << std::fixed << gps_time << ' ' << position << ' ' << velocity << '\n';
    }

    file.close();

    return true_states;
}

std::vector<RawMeasurementGroupped> DataParser::load_grace_gps_data(const std::string& filename) {
    std::vector<RawMeasurementGroupped> raw_mgs;
    std::vector<RawMeasurement> raw_ms(32);
    int gps_time_current = 0;
    
    std::ifstream file;
    file.open("../data/grace/gps/" + filename, std::ios::binary);

    std::string line;
    unsigned skiprows = 24;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    char record[74];
    while (!file.eof()) {
        file.read(record, 4); 
        int gps_time = bytes_to_T_endian(record, int());
        
        if (static_cast<int>(record[0]) == 0) break;

        file.read(record, 5);

        file.read(record, 1);
        unsigned prn_id = static_cast<unsigned>(record[0]);

        file.read(record, 3);
        
        file.read(record, 1);
        unsigned char qualflg = static_cast<unsigned char>(record[0]);

        file.read(record, 8);

        file.read(record, 8);
        double L1_range = bytes_to_T_endian(record, double());
        file.read(record, 8);
        double L2_range = bytes_to_T_endian(record, double());

        file.read(record, 8);

        file.read(record, 8);
        double L1_phase = bytes_to_T_endian(record, double());
        file.read(record, 8);
        double L2_phase = bytes_to_T_endian(record, double());

        file.read(record, 2);

        file.read(record, 2);
        double L1_SNR = bytes_to_T_endian(record, double());
        file.read(record, 2);
        double L2_SNR = bytes_to_T_endian(record, double());

        file.read(record, 6);

        RawMeasurement raw_m = {true,
                                gps_time, prn_id,
                                L1_range, L2_range,
                                L1_phase, L2_phase,
                                L1_SNR, L2_SNR,
                                qualflg};

        if (gps_time_current == 0) gps_time_current = gps_time;

        if (gps_time != gps_time_current) {
            raw_mgs.push_back({gps_time_current, raw_ms});
            
            raw_ms.clear();
            raw_ms.resize(32);

            gps_time_current = gps_time;
        }

        size_t prn_index = prn_id - 1;
        raw_ms[prn_index] = raw_m;
    }

    file.close();

    return raw_mgs;
}

std::vector<AccelerationMeasurement> DataParser::load_grace_acc_data(const std::string& filename) {
    std::vector<AccelerationMeasurement> acc_ms;

    std::ifstream file;
    file.open("../data/grace/acc/" + filename, std::ios::binary);

    std::string line;
    unsigned skiprows = 24;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    char record[78];
    while (!file.eof()) {
        file.read(record, 4);
        int gps_time = bytes_to_T_endian(record, int());

        if (static_cast<int>(record[0]) == 0) break;

        file.read(record, 1);

        file.read(record, 8);
        double lin_accl_x = bytes_to_T_endian(record, double());
        file.read(record, 8);
        double lin_accl_y = bytes_to_T_endian(record, double());
        file.read(record, 8);
        double lin_accl_z = bytes_to_T_endian(record, double());

        file.read(record, 8);
        double ang_accl_x = bytes_to_T_endian(record, double());
        file.read(record, 8);
        double ang_accl_y = bytes_to_T_endian(record, double());
        file.read(record, 8);
        double ang_accl_z = bytes_to_T_endian(record, double());

        file.read(record, 8);
        double acl_x_res = bytes_to_T_endian(record, double());
        file.read(record, 8);
        double acl_y_res = bytes_to_T_endian(record, double());
        file.read(record, 8);
        double acl_z_res = bytes_to_T_endian(record, double());

        file.read(record, 1);

        std::vector<double> lin_acc = {lin_accl_x, lin_accl_y, lin_accl_z};
        std::vector<double> ang_acc = {ang_accl_x, ang_accl_y, ang_accl_z};
        acc_ms.push_back({gps_time, lin_acc, ang_acc});
    }

    file.close();

    return acc_ms;
}

std::vector<State> DataParser::load_swarm_nav_data(const std::string& filename) {
    std::vector<State> true_states;

    std::ifstream file;
    file.open("../data/swarm/nav/" + filename, std::ios::in);

    std::string line;

    unsigned skiprows = 22;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    int count = 0;
    while (!file.eof()) {
        std::getline(file, line);

        if (line[0] == 'E') break;

        int year = std::stoi(line.substr(3, 4));
        int month = std::stoi(line.substr(8, 2));
        int day = std::stoi(line.substr(11, 2));
        int hours = std::stoi(line.substr(14, 2));
        int minutes = std::stoi(line.substr(17, 2));
        int seconds = std::stoi(line.substr(20, 2));

        std::getline(file, line);

        double xpos = std::stod(line.substr(5, 13));
        double ypos = std::stod(line.substr(19, 13));
        double zpos = std::stod(line.substr(33, 13));

        std::getline(file, line);

        double xvel = std::stod(line.substr(5, 13));
        double yvel = std::stod(line.substr(19, 13));
        double zvel = std::stod(line.substr(33, 13));

        int seconds_since_midnight = seconds + minutes * 60 + hours * 3600;

        std::vector<double> position = {xpos, ypos, zpos};
        std::vector<double> velocity = {xvel, yvel, zvel};
        true_states.push_back({seconds_since_midnight, position, velocity}); // Перевести в gps_time?
    }

    file.close();

    return true_states;
}

std::vector<RawMeasurementGroupped> DataParser::load_swarm_gps_data(const std::string& filename) {
    std::vector<RawMeasurementGroupped> raw_mgs;
    std::vector<RawMeasurement> raw_ms(32);
    int gps_time_current = 0;
    
    std::ifstream file;
    file.open("../data/swarm/gps/" + filename, std::ios::binary);

    std::string line;

    unsigned skiprows = 17;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {

        int year = std::stoi(line.substr(2, 4));
        int month = std::stoi(line.substr(7, 2));
        int day = std::stoi(line.substr(10, 2));
        int hours = std::stoi(line.substr(13, 2));
        int minutes = std::stoi(line.substr(16, 2));
        int seconds = std::stoi(line.substr(19, 2));

        int seconds_since_midnight = seconds + minutes * 60 + hours * 3600;

        int sat_count = std::stoi(line.substr(33, 4));

        for (size_t i = 0; i < sat_count; i++) {
            std::getline(file, line);

            unsigned prn_id = std::stoi(line.substr(1, 2));
            double C1C = std::stod(line.substr(3, 14));
            double L1C = std::stod(line.substr(19, 14));
            double S1C = std::stod(line.substr(35, 14));
            double C1W = std::stod(line.substr(51, 14));
            double S1W = std::stod(line.substr(67, 14));
            double C2W = std::stod(line.substr(83, 14));
            double L2W = std::stod(line.substr(99, 14));
            double S2W = std::stod(line.substr(115, 14));

            RawMeasurement raw_m = {true,
                                    seconds_since_midnight, prn_id,
                                    C1C, C2W,
                                    L1C, L2W,
                                    S1C, S2W, 
                                    0};

            size_t prn_index = prn_id - 1;
            raw_ms[prn_index] = raw_m;
        }

        raw_mgs.push_back({seconds_since_midnight, raw_ms});
            
        raw_ms.clear();
        raw_ms.resize(32);
    }

    return raw_mgs;
}

std::vector<std::vector<Ephemeris>> DataParser::load_brdc_data(const std::string& filename) {
    std::vector<std::vector<Ephemeris>> ephs(32);
    Ephemeris eph;

    std::ifstream file;
    file.open("../data/gps/brdc/" + filename, std::ios::in);

    std::string line;

    unsigned skiprows = 8;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {
        eph.prn_id = std::stoi(line.substr(0, 2));
        eph.a_f0 = std::stod(line.substr(22, 19));
        eph.a_f1 = std::stod(line.substr(41, 19));
        eph.a_f2 = std::stod(line.substr(60, 19));

        std::getline(file, line);
 
        eph.C_rs = std::stod(line.substr(22, 19));
        eph.delta_n = std::stod(line.substr(41, 19));
        eph.M_0 = std::stod(line.substr(60, 19));

        std::getline(file, line);

        eph.C_uc = std::stod(line.substr(3, 19));
        eph.e = std::stod(line.substr(22, 19));
        eph.C_us = std::stod(line.substr(41, 19));
        eph.A_sqrt = std::stod(line.substr(60, 19));

        std::getline(file, line);

        eph.t_oe = static_cast<int>(std::stod(line.substr(3, 19)));
        eph.C_ic = std::stod(line.substr(22, 19));
        eph.Omega_0 = std::stod(line.substr(41, 19));
        eph.C_is = std::stod(line.substr(60, 19));
    
        std::getline(file, line);

        eph.i_0 = std::stod(line.substr(3, 19));
        eph.C_rc = std::stod(line.substr(22, 19));
        eph.omega = std::stod(line.substr(41, 19));
        eph.Omega_dot = std::stod(line.substr(60, 19));

        std::getline(file, line);

        eph.IDOT = std::stod(line.substr(3, 19));

        std::getline(file, line);
        std::getline(file, line);

        unsigned prn_index = eph.prn_id - 1;
        ephs[prn_index].push_back(eph);
    }

    file.close();

    return ephs;
}

