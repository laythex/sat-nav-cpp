#pragma once

class Config {
public:
    Config() = delete;

    inline static double GDOP0;
    inline static double CN0_min_threshold;
    inline static double CN0_max_threshold;
    inline static unsigned fadein_threshold;
    inline static double hatch_constant;
    inline static double solution_tolerance;
};