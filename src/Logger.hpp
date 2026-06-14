#pragma once

#include <fstream>
#include <string_view>

#include "Structures.hpp"

class Logger {
public:
    Logger(std::string_view filename);

    void lnbr(char end = '\n');
    void log(std::string_view arg, char end = '\n');

    template <typename T> void log(const T& arg, char end = '\n') {
        if constexpr (requires { to_string(arg); }) {
            file << to_string(arg) << end;
        } else {
            file << arg << end;
        }
    }

    template <typename T> void logv(std::string_view desc, const T& val, char end = '\n') {
        file << desc << ":\t";
        if constexpr (requires { to_string(val); }) {
            file << to_string(val) << end;
        } else {
            file << val << end;
        }
    }

private:
    std::ofstream file;
};
