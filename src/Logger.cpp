#include "Logger.hpp"

#include <filesystem>
#include <format>
#include <iostream>
#include <stdexcept>

Logger::Logger(std::string_view filename) {
    std::filesystem::path dir("../logs");
    std::filesystem::create_directories(dir);

    std::filesystem::path filepath = dir / std::format("{}.txt", filename);
    file.open(filepath, std::ofstream::trunc);
}

void Logger::lnbr(char end) { file << end; }

void Logger::log(std::string_view arg, char end) { file << arg << end; }
