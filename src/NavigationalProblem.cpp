#include "NavigationalProblem.hpp"

NavigationalProblem::NavigationalProblem(std::string gnv_filename, std::string gps_filename) {
    load_gnv_data(gnv_filename);
    load_gps_data(gps_filename);
}

void NavigationalProblem::load_gnv_data(std::string gnv_filename) {
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

void NavigationalProblem::load_gps_data(std::string gps_filename) {
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
        std::pair key = {prn_id, gps_time};
        pr[key] = {L1_range, L2_range};
    }

    gps_file.close();
}

void NavigationalProblem::solve(unsigned ti, unsigned tf) {

    std::set<unsigned>::iterator iti = ts.lower_bound(ti);
    std::set<unsigned>::iterator itf = ts.upper_bound(tf);
    for (auto it = iti; it != itf; it++) {
        unsigned t = *it;

        get_pr_and_gps_pos(t);
    }
}