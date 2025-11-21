#include "DataParser.hpp"

template <typename T>
T bytes_to_T_endian(const char* bytes, T t) {
    size_t n = sizeof(T);
    char res_bytes[n];

    for (unsigned i = 0; i < n; i++) {
        res_bytes[n - i - 1] = bytes[i];
    }

    return *reinterpret_cast<T*>(res_bytes);
}

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
    file.open("../data/gps/" + filename, std::ios::in);

    std::string line;

    unsigned skiprows = 196;
    for (unsigned i = 0; i < skiprows; ++i) {
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

        unsigned short L1_SNR = std::stoi(line);
        
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        unsigned short L2_SNR = std::stoi(line);

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

std::vector<State> DataParser::load_grace_gnv_data(const std::string& filename) {
    std::vector<State> true_states;

    std::ifstream file;
    file.open("../data/gnv/" + filename, std::ios::binary);

    std::string line;
    unsigned skiprows = 26;
    for (unsigned i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

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
    }

    file.close();

    return true_states;
}

std::vector<RawMeasurementGroupped> DataParser::load_grace_gps_data(const std::string& filename) {
    std::vector<RawMeasurementGroupped> raw_mgs;
    std::vector<RawMeasurement> raw_ms(32);
    int gps_time_current = 0;
    
    std::ifstream file;
    file.open("../data/gps/" + filename, std::ios::binary);

    std::string line;
    unsigned skiprows = 24;
    for (unsigned i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    char record[74];
    while (!file.eof()) {
        file.read(record, 4); 
        int gps_time = bytes_to_T_endian(record, int()); // сделать микросекунды?
        
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
        unsigned short L1_SNR = bytes_to_T_endian(record, (unsigned short)(0));
        file.read(record, 2);
        unsigned short L2_SNR = bytes_to_T_endian(record, (unsigned short)(0));

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

        unsigned prn_index = prn_id - 1;
        raw_ms[prn_index] = raw_m;
    }

    file.close();

    return raw_mgs;
}

std::vector<std::vector<Ephemeris>> DataParser::load_brdc_data(const std::string& filename) {
    std::vector<std::vector<Ephemeris>> ephs(32);
    Ephemeris eph;

    std::ifstream file;
    file.open("../data/brdc/" + filename, std::ios::in);

    std::string line;

    unsigned skiprows = 8;
    for (unsigned i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {
        eph.prn_id = std::stoi(line.substr(0, 2));
        eph.a_f0 = std::stod(line.substr(22, 41));
        eph.a_f1 = std::stod(line.substr(41, 60));
        eph.a_f2 = std::stod(line.substr(60, 79));

        std::getline(file, line);

        eph.C_rs = std::stod(line.substr(22, 41));
        eph.delta_n = std::stod(line.substr(41, 60));
        eph.M_0 = std::stod(line.substr(60, 79));

        std::getline(file, line);

        eph.C_uc = std::stod(line.substr(3, 22));
        eph.e = std::stod(line.substr(22, 41));
        eph.C_us = std::stod(line.substr(41, 60));
        eph.A_sqrt = std::stod(line.substr(60, 79));

        std::getline(file, line);

        eph.t_oe = static_cast<int>(std::stod(line.substr(3, 22)));
        eph.C_ic = std::stod(line.substr(22, 41));
        eph.Omega_0 = std::stod(line.substr(41, 60));
        eph.C_is = std::stod(line.substr(60, 79));
    
        std::getline(file, line);

        eph.i_0 = std::stod(line.substr(3, 22));
        eph.C_rc = std::stod(line.substr(22, 41));
        eph.omega = std::stod(line.substr(41, 60));
        eph.Omega_dot = std::stod(line.substr(60, 79));

        std::getline(file, line);

        eph.IDOT = std::stod(line.substr(3, 22));

        std::getline(file, line);
        std::getline(file, line);

        unsigned prn_index = eph.prn_id - 1;
        ephs[prn_index].push_back(eph);
    }

    file.close();

    return ephs;
}
