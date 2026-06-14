#include "DataParser.hpp"

std::vector<State> DataParser::load_grace_fo_gnv_data(const Date& date, GRACE sat) {
    std::vector<State> true_states;

    std::string filename = grace_gnv_dir + "GNV1B_" + date_to_grace_string(date) + "_" + to_char(sat) + "_04.txt";

    std::ifstream file;
    file.open(filename, std::ios::in);

    std::string line;

    unsigned skiprows = 148;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {
        std::stringstream stream(line);

        std::getline(stream, line, ' ');
        unsigned grace_time = static_cast<unsigned>(std::stoul(line));
        unsigned gps_time = grace_to_gps_time(grace_time);

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

        Eigen::Vector3d position = {xpos, ypos, zpos};
        Eigen::Vector3d velocity = {xvel, yvel, zvel};
        true_states.push_back({gps_time, position, velocity});
    }

    file.close();

    return true_states;
}

std::vector<RawMeasurementGroupped> DataParser::load_grace_fo_gps_data(const Date& date, GRACE sat) {
    std::vector<RawMeasurementGroupped> raw_mgs;
    std::vector<RawMeasurement> raw_ms(32);
    unsigned gps_time_current = 0;

    std::string filename = grace_gps_dir + "GPS1B_" + date_to_grace_string(date) + "_" + to_char(sat) + "_04.txt";

    std::ifstream file;
    file.open(filename, std::ios::in);

    std::string line;

    unsigned skiprows = 196;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {
        std::stringstream stream(line);

        std::getline(stream, line, ' ');

        unsigned grace_time = static_cast<unsigned>(std::stoul(line));
        unsigned gps_time = grace_to_gps_time(grace_time);

        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        unsigned prn_id = static_cast<unsigned>(std::stoul(line));

        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        unsigned qualflg = static_cast<unsigned>(std::stoul(line));

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

        double L1_CN0 = convert_SNR_to_CN0(std::stod(line));

        std::getline(stream, line, ' ');
        std::getline(stream, line, ' ');

        double L2_CN0 = convert_SNR_to_CN0(std::stod(line));

        RawMeasurement raw_m = {true,     gps_time, prn_id, L1_range, L2_range,
                                L1_phase, L2_phase, L1_CN0, L2_CN0,   qualflg};

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

std::vector<State> DataParser::load_grace_gnv_data(const Date& date, GRACE sat) {
    std::vector<State> true_states;

    std::string filename = grace_gnv_dir + "GNV1B_" + date_to_grace_string(date) + "_" + to_char(sat) + "_02.dat";

    std::ifstream file;
    file.open(filename, std::ios::binary);

    std::string line;
    unsigned skiprows = 26;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    char record[103];
    while (!file.eof()) {
        file.read(record, 4);
        unsigned grace_time = bytes_to_T_endian(record, unsigned());
        unsigned gps_time = grace_to_gps_time(grace_time);

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

        Eigen::Vector3d position = {xpos, ypos, zpos};
        Eigen::Vector3d velocity = {xvel, yvel, zvel};
        true_states.push_back({gps_time, position, velocity});
    }

    file.close();

    return true_states;
}

std::vector<RawMeasurementGroupped> DataParser::load_grace_gps_data(const Date& date, GRACE sat) {
    std::vector<RawMeasurementGroupped> raw_mgs;
    std::vector<RawMeasurement> raw_ms(32);
    unsigned gps_time_current = 0;

    std::string filename = grace_gps_dir + "GPS1B_" + date_to_grace_string(date) + "_" + to_char(sat) + "_02.dat";

    std::ifstream file;
    file.open(filename, std::ios::binary);

    std::string line;
    unsigned skiprows = 24;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    char record[74];
    while (!file.eof()) {
        file.read(record, 4);
        unsigned grace_time = bytes_to_T_endian(record, unsigned());
        unsigned gps_time = grace_to_gps_time(grace_time);

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
        double L1_CN0 = convert_SNR_to_CN0(bytes_to_T_endian(record, (unsigned short)(0)));
        file.read(record, 2);
        double L2_CN0 = convert_SNR_to_CN0(bytes_to_T_endian(record, (unsigned short)(0)));

        file.read(record, 6);

        RawMeasurement raw_m = {true,     gps_time, prn_id, L1_range, L2_range,
                                L1_phase, L2_phase, L1_CN0, L2_CN0,   qualflg};

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

std::vector<AccelerationMeasurement> DataParser::load_grace_acc_data(const Date& date, GRACE sat) {
    std::vector<AccelerationMeasurement> acc_ms;

    std::string filename = grace_acc_dir + "GNV1B_" + date_to_grace_string(date) + "_" + to_char(sat) + "_02.dat";

    std::ifstream file;
    file.open(filename, std::ios::binary);

    std::string line;
    unsigned skiprows = 24;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    char record[78];
    while (!file.eof()) {
        file.read(record, 4);
        unsigned grace_time = bytes_to_T_endian(record, unsigned());
        unsigned gps_time = grace_to_gps_time(grace_time);

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
        file.read(record, 8);
        file.read(record, 8);

        file.read(record, 1);

        std::vector<double> lin_acc = {lin_accl_x, lin_accl_y, lin_accl_z};
        std::vector<double> ang_acc = {ang_accl_x, ang_accl_y, ang_accl_z};
        acc_ms.push_back({gps_time, lin_acc, ang_acc});
    }

    file.close();

    return acc_ms;
}

std::vector<State> DataParser::load_swarm_nav_data(const Date& date, SWARM sat) {
    std::vector<State> true_states;

    std::string filename = swarm_nav_dir + "SW_OPER_GPS" + to_char(sat) + "NAV_1B_" + date_to_swarm_string(date) +
                           "T000000_" + date_to_swarm_string(date) + "T235959_0602" + ".sp3";
    std::ifstream file;
    file.open(filename, std::ios::in);

    std::string line;

    unsigned skiprows = 22;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (!file.eof()) {
        std::getline(file, line);

        if (line[0] == 'E') break;

        unsigned year = static_cast<unsigned>(std::stoul(line.substr(3, 4)));
        unsigned month = static_cast<unsigned>(std::stoul(line.substr(8, 2)));
        unsigned day = static_cast<unsigned>(std::stoul(line.substr(11, 2)));
        unsigned hours = static_cast<unsigned>(std::stoul(line.substr(14, 2)));
        unsigned minutes = static_cast<unsigned>(std::stoul(line.substr(17, 2)));
        unsigned seconds = static_cast<unsigned>(std::stoul(line.substr(20, 2)));

        unsigned gps_time = date_to_gps_time(year, month, day, hours, minutes, seconds);

        std::getline(file, line);

        double xpos = std::stod(line.substr(5, 13)) * 1e3;
        double ypos = std::stod(line.substr(19, 13)) * 1e3;
        double zpos = std::stod(line.substr(33, 13)) * 1e3;

        std::getline(file, line);

        double xvel = std::stod(line.substr(5, 13));
        double yvel = std::stod(line.substr(19, 13));
        double zvel = std::stod(line.substr(33, 13));

        Eigen::Vector3d position = {xpos, ypos, zpos};
        Eigen::Vector3d velocity = {xvel, yvel, zvel};
        true_states.push_back({gps_time, position, velocity});
    }

    file.close();

    return true_states;
}

std::vector<RawMeasurementGroupped> DataParser::load_swarm_gps_data(const Date& date, SWARM sat) {
    std::vector<RawMeasurementGroupped> raw_mgs;
    std::vector<RawMeasurement> raw_ms(32);

    std::string filename = swarm_gps_dir + "SW_OPER_GPS" + to_char(sat) + "_RO_1B_" + date_to_swarm_string(date) +
                           "T000000_" + date_to_swarm_string(date) + "T235959_0602" + ".rnx";

    std::ifstream file;
    file.open(filename, std::ios::binary);

    std::string line;

    unsigned skiprows = 17;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {

        unsigned year = static_cast<unsigned>(std::stoul(line.substr(2, 4)));
        unsigned month = static_cast<unsigned>(std::stoul(line.substr(7, 2)));
        unsigned day = static_cast<unsigned>(std::stoul(line.substr(10, 2)));
        unsigned hours = static_cast<unsigned>(std::stoul(line.substr(13, 2)));
        unsigned minutes = static_cast<unsigned>(std::stoul(line.substr(16, 2)));
        unsigned seconds = static_cast<unsigned>(std::stoul(line.substr(19, 2)));

        unsigned gps_time = date_to_gps_time(year, month, day, hours, minutes, seconds);

        size_t sat_count = std::stoul(line.substr(33, 4));

        for (size_t i = 0; i < sat_count; ++i) {
            std::getline(file, line);

            unsigned prn_id = static_cast<unsigned>(std::stoul(line.substr(1, 2)));

            double C1C = std::stod(line.substr(3, 14));  // L1 C/A PR
            double L1C = std::stod(line.substr(19, 14)); // L1 C/A CP
            double S1C = std::stod(line.substr(35, 14)); // L1 C/A CN0
            // double C1W = std::stod(line.substr(51, 14));    // L1 Z   PR
            // double S1W = std::stod(line.substr(67, 14));    // L1 Z   CN0
            double C2W = std::stod(line.substr(83, 14));  // L2 Z   PR
            double L2W = std::stod(line.substr(99, 14));  // L2 Z   CP
            double S2W = std::stod(line.substr(115, 14)); // L2 Z   CN0
            // std::cout << std::fixed << seconds << ' ' << prn_id << ' ' << C1C << ' ' << C2W << ' ' << L2W << ' ' <<
            // S2W << std::endl;
            RawMeasurement raw_m = {true, gps_time, prn_id, C1C,
                                    C2W,  L1C,      L2W, // Отсутствует L1W, поэтому hatch-фильтр не используем
                                    S1C,  S2W,      0};

            size_t prn_index = prn_id - 1;
            raw_ms[prn_index] = raw_m;
        }

        raw_mgs.push_back({gps_time, raw_ms});

        raw_ms.clear();
        raw_ms.resize(32);
    }

    return raw_mgs;
}

// Заменять D на E
std::vector<std::vector<Ephemeris>> DataParser::load_brdc_data(const Date& date) {
    std::vector<std::vector<Ephemeris>> ephs(32);
    Ephemeris eph;

    std::string filename = gps_brdc_dir + date_to_gps_string(date);

    std::ifstream file;
    file.open(filename, std::ios::in);

    std::string line;

    unsigned skiprows = 8;
    for (size_t i = 0; i < skiprows; ++i) {
        std::getline(file, line);
    }

    while (std::getline(file, line)) {
        eph.prn_id = static_cast<unsigned>(std::stoul(line.substr(0, 2)));
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

unsigned DataParser::grace_to_gps_time(unsigned grace_time) { return grace_time + 630763200; }

unsigned DataParser::date_to_gps_time(unsigned year, unsigned month, unsigned day, unsigned hours, unsigned minutes,
                                      unsigned seconds) {
    unsigned day_of_the_year = get_day_of_the_year({year, month, day});
    unsigned seconds_since_midnight = seconds + minutes * 60 + hours * 3600;

    return (year - 1980) * 365 * day_seconds_count + (day_of_the_year + 6) * day_seconds_count + seconds_since_midnight;
}

std::string DataParser::date_to_grace_string(const Date& date) {
    return std::format("{:04}-{:02}-{:02}", date.year, date.month, date.day);
}

std::string DataParser::date_to_swarm_string(const Date& date) {
    return std::format("{:04}{:02}{:02}", date.year, date.month, date.day);
}

std::string DataParser::date_to_gps_string(const Date& date) {
    unsigned day_of_the_year = get_day_of_the_year(date);
    return std::format("brdc{:03}0.{:02}n", day_of_the_year, date.year % 100);
}

bool DataParser::is_leap_year(unsigned year) { return (year % 4 == 0) && (year % 100 != 0 || year % 400 == 0); }

unsigned DataParser::get_day_of_the_year(const Date& date) {
    unsigned day_of_the_year = 0;

    for (size_t i = 0; i < date.month - 1; ++i) {
        day_of_the_year += month_day_count[i];

        if (i == 1 && is_leap_year(date.year)) {
            ++day_of_the_year;
        }
    }

    return day_of_the_year + date.day;
}

double DataParser::convert_SNR_to_CN0(double SNR) { return 20 * log10(SNR) - 3.01029995664; }
