#include "Logger.hpp"

#include <iostream>

Logger::Logger(const std::string& filename) {
    file.open(std::format("../logs/{}.txt", filename), std::ofstream::trunc);
}

Logger::~Logger() {
    file.close();
}

void Logger::lnbr(char end) {
    file << end;
}

void Logger::log(const std::string& arg, char end) {
    file << arg << end;
}

void Logger::log(double arg, char end) {
    file << std::to_string(arg) << end;
}

void Logger::log(size_t arg, char end) {
    file << std::to_string(arg) << end;
}

void Logger::logv(const std::string& desc, double val, char end) {
    file << desc << ":\t" << std::to_string(val) << end;
}

void Logger::logv(const std::string& desc, size_t val, char end) {
    file << desc << ":\t" << std::to_string(val) << end;
}

void Logger::logv(const std::string& desc, int val, char end) {
    file << desc << ":\t" << std::to_string(val) << end;
}

void Logger::logv(const std::string& desc, unsigned val, char end) {
    file << desc << ":\t" << std::to_string(val) << end;
}