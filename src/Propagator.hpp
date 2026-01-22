#pragma once

#include <cmath>
#include <vector>

#include "LinAlg.hpp"
#include "Structures.hpp"

namespace Propagator {

State propagate_rk4(const State& state, double timespan);

std::vector<double> rhs(const std::vector<double>& arg);

};