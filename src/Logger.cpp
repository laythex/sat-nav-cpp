#include "Logger.hpp"

Logger::Logger() {
    file.open("../logs/log.txt", std::fstream::out);
}

Logger::~Logger() {
    file.close();
}

void Logger::log(const std::string& line) {
    file << line << std::endl;
}