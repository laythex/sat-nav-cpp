#pragma once

#include <string>
#include <format>
#include <fstream>

class Logger {
public:
    Logger(const std::string& filename);
    ~Logger();
    
    void lnbr(char end = '\n');

    void log(const std::string& arg, char end = '\n');
    void log(double arg, char end = '\n');
    void log(size_t arg, char end = '\n');

    void logv(const std::string& desc, double val, char end = '\n');
    void logv(const std::string& desc, size_t val, char end = '\n');
    void logv(const std::string& desc, int val, char end = '\n');
    void logv(const std::string& desc, unsigned val, char end = '\n');

private:
    std::ofstream file;
};