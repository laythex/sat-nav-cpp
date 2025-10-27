#pragma once

class GPSHandler {

public:
    GPSHandler(string rnx_data_filename);

    std::vector<double> interp_state(unsigned prn_id, unsigned gps_time);

    static unsigned convert_id(unsigned prn_id);

private:
    void load_rnx_data(string rnx_data_filename);

    double mu = 3.986005e14;
    double Omega_e_dot = 7.2921151467e-5;
    double F = -4.442807633e-10

    std::map<std::pair<unsigned, unsigned>, double> map_M_0;
    std::map<std::pair<unsigned, unsigned>, double> map_delta_n;
    std::map<std::pair<unsigned, unsigned>, double> map_e;
    std::map<std::pair<unsigned, unsigned>, double> map_A_sqrt;
    std::map<std::pair<unsigned, unsigned>, double> map_Omega_0;
    std::map<std::pair<unsigned, unsigned>, double> map_i_0;
    std::map<std::pair<unsigned, unsigned>, double> map_omega;
    std::map<std::pair<unsigned, unsigned>, double> map_Omega_dot;
    std::map<std::pair<unsigned, unsigned>, double> map_IDOT;
    std::map<std::pair<unsigned, unsigned>, double> map_C_uc;
    std::map<std::pair<unsigned, unsigned>, double> map_C_us;
    std::map<std::pair<unsigned, unsigned>, double> map_C_rc;
    std::map<std::pair<unsigned, unsigned>, double> map_C_rs;
    std::map<std::pair<unsigned, unsigned>, double> map_C_ic;
    std::map<std::pair<unsigned, unsigned>, double> map_C_is;
    std::map<std::pair<unsigned, unsigned>, double> map_t_oe;

    std::map<std::pair<unsigned, unsigned>, double> map_a_f0;
    std::map<std::pair<unsigned, unsigned>, double> map_a_f1;
    std::map<std::pair<unsigned, unsigned>, double> map_a_f2;
}