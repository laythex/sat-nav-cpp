#pragma once

#include <cmath>
#include <vector>

#include "LinAlg.hpp"
#include "Structures.hpp"

namespace Propagator { // переделать в статик класс?

State propagate_rk4(const State& state, double timespan);

State propagate_rel_rk4(const State& init_state, double timespan, const std::vector<double>& init_pos_pas, const std::vector<double>& fin_pos_pas);

std::vector<double> rhs(const std::vector<double>& arg);

std::vector<double> rhs_rel(const std::vector<double>& arg, const std::vector<double>& pos_pas);

};