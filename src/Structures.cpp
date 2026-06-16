#include "Structures.hpp"

bool StandaloneData::is_fully_solved() const {
    bool pas_solved = pas.status == SolutionStatus::OK;
    bool act_solved = act.status == SolutionStatus::OK;

    return pas_solved && act_solved;
}

bool StandaloneData::is_fully_present_at(size_t prn_id) const {
    size_t prn_index = prn_id - 1;

    bool gps_pas_present = pas.gps_mask[prn_index];
    bool gps_act_present = act.gps_mask[prn_index];

    return gps_pas_present && gps_act_present;
}