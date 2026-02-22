#include "Propagator.hpp"

State Propagator::propagate_rk4(const State& init_state, double timespan) {
    std::vector<double> arg = {init_state.position[0], init_state.position[1], init_state.position[2],
                               init_state.velocity[0], init_state.velocity[1], init_state.velocity[2]};

    double tau = timespan / 1e4; // Определить оптимальный шаг
    double time_current = 0;

    std::vector<double> k1(6), k2(6), k3(6), k4(6);

    while (true) {
        std::vector<double> k1 = rhs(arg);
        std::vector<double> k2 = rhs(arg + k1 * (tau / 2));
        std::vector<double> k3 = rhs(arg + k2 * (tau / 2));
        std::vector<double> k4 = rhs(arg + k3 * tau);

        arg = arg + (k1 + k2 * 2 + k3 * 2 + k4) * (tau / 6);

        time_current += tau;
        if (time_current >= timespan) break;
    }

    return {init_state.time, {arg[0], arg[1], arg[2]}, {arg[3], arg[4], arg[5]}};
}

State Propagator::propagate_rel_rk4(const State& init_state, double timespan, const std::vector<double>& init_pos_pas, const std::vector<double>& fin_pos_pas) {
    std::vector<double> arg = {init_state.position[0], init_state.position[1], init_state.position[2],
                               init_state.velocity[0], init_state.velocity[1], init_state.velocity[2]};

    double tau = timespan / 1e2; // Определить оптимальный шаг
    double time_current = 0;

    std::vector<double> k1(6), k2(6), k3(6), k4(6);

    while (true) {
        std::vector<double> k1 = rhs_rel(arg, init_pos_pas);
        std::vector<double> k2 = rhs_rel(arg + k1 * (tau / 2), (init_pos_pas + fin_pos_pas) / 2);
        std::vector<double> k3 = rhs_rel(arg + k2 * (tau / 2), (init_pos_pas + fin_pos_pas) / 2);
        std::vector<double> k4 = rhs_rel(arg + k3 * tau, fin_pos_pas);

        arg = arg + (k1 + k2 * 2 + k3 * 2 + k4) * (tau / 6);

        time_current += tau;
        if (time_current >= timespan) break;
    }

    return {init_state.time, {arg[0], arg[1], arg[2]}, {arg[3], arg[4], arg[5]}};
}

std::vector<double> Propagator::rhs(const std::vector<double>& arg) {

    double r = sqrt(arg[0] * arg[0] + arg[1] * arg[1] + arg[2] * arg[2]);
    double k = -3.986004418e14 / (r * r * r);

    return {arg[3], arg[4], arg[5], k * arg[0], k * arg[1], k * arg[2]};
}

std::vector<double> Propagator::rhs_rel(const std::vector<double>& arg, const std::vector<double>& pos_pas) {
 
    std::vector<double> pos = {arg[0], arg[1], arg[2]};
    std::vector<double> vel = {arg[3], arg[4], arg[5]};

    double r_pas = abs(pos_pas); 
    double k = 3 * -3.986004418e14 / (r_pas * r_pas * r_pas * r_pas * r_pas);
    std::vector<double> acc = col_by_row(pos_pas, pos_pas) * pos * k;

    return {vel[0], vel[1], vel[2], acc[0], acc[1], acc[2]};
}
