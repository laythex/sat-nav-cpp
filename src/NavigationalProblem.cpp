#include "NavigationalProblem.hpp"

NavigationalProblem::NavigationalProblem(std::string gnv_filename, std::string gps_filename, const GPSHandler& handler) : handler(handler) {
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
        pr[{prn_id, gps_time}] = {L1_range, L2_range};
    }

    gps_file.close();
}

void NavigationalProblem::solve(unsigned ti, unsigned tf) {

    auto it_i = ts.lower_bound(ti);
    auto it_f = ts.upper_bound(tf);
    for (auto it = it_i; it != it_f; it++) {
        unsigned t = *it;

        std::vector<double> pseudoranges;
        std::map<unsigned, std::vector<double>> gps_positions;

        for (unsigned prn_id = 1; prn_id <= 32; prn_id++) {
            auto it = pr.find({prn_id, t});
            if (it == pr.end()) continue;

            std::vector<double> measurments = it->second;
            double L1_range = measurments[0];
            double L2_range = measurments[1];

            double gps_time = t;
            double pseudorange = C1 * L1_range + C2 * L2_range;

            //     Эффект                     t       t * v    t * c   
            // (1) Задержка распространения - 80 мс - 240 км -  ---   
            // (2) Ошибка часов НКА         - 3 мс  - 10 м   - 1000 км 
            // (3) Релятивизм               - 1 мкс - 3 мм   - 300 м  

            // При расчете ошибки часов учитываем только (1)
            // При расчете эфемерид учитываем (1), (2) 
            // При корректировке псевдодальности учитываем (1), (2), (3)

            double propagation_delay = pseudorange / c;
            gps_time -= propagation_delay;

            double clock_error = handler.get_clock_error(prn_id, gps_time);
            gps_time -= clock_error;

            std::vector<double> state = handler.get_state(prn_id, gps_time);

            double relativity_error = state[6];
            gps_time -= relativity_error;

            pseudorange += (clock_error + relativity_error) * c;

            pseudoranges[prn_id] = pseudorange;
            gps_positions[prn_id] = std::vector<double>(state.begin(), state.begin() + 3);
        }

        iterative(pseudoranges, gps_positions);
    }
}

void NavigationalProblem::iterative(const std::vector<double>& PR, const std::map<unsigned, std::vector<double>>& X) {
    
    std::vector<double> x0 = {0, 0, 0};

    // while (true) {


    // }
}
