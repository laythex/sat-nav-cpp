#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "DataParser.hpp"
#include "GPSHandler.hpp"
#include "Logger.hpp"
#include "SatNav.hpp"
#include "SatNavUtils.hpp"
#include "Structures.hpp"

class SatNav;

class SatNavRel {

public:
    SatNavRel(SatNav& passive, SatNav& active);
    SatNavRel(const SatNavRel& sn);

    void solve(unsigned ti = 0, unsigned tf = 0);

    const std::vector<State>& get_true_states() const;
    const std::vector<SolutionStateRel>& get_solution_states() const;

    std::vector<State>::const_iterator get_true_state_iterator(unsigned time) const;

    // Перенести в приват
    SatNav& pas;
    SatNav& act;

private:
    enum class SatType { PASSIVE, ACTIVE };

    void solve_standalones(unsigned ti = 0, unsigned tf = 0);
    void form_true_states();
    void form_standalone_data();

    SolutionStateRel find_solution_ls(const StandaloneData& data);

    std::vector<State> true_states;
    std::vector<StandaloneData> standalone_datas;
    std::vector<SolutionStateRel> solution_states;

    Logger logger;
};
