#pragma once

class Config {
public:
    Config() = delete;

    inline static unsigned ti;
    inline static unsigned tf;

    inline static double GDOP0;
    inline static double CN0_min_threshold;
    inline static double CN0_max_threshold;
    inline static unsigned fadein_threshold;
    inline static double hatch_constant;
    inline static double solution_tolerance;

    inline static bool use_kf;
    inline static double residual_threshold;
    inline static double residual_rel_threshold;
    inline static double lambda_r;
    inline static double lambda_v;
    inline static double lambda_model;
    inline static double trace_threshold;
    inline static double use_rk4;
};