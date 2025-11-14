#pragma once

#include <string>

#include <fstream>

class Logger {
public:
    Logger();
    ~Logger();
    
    void log(const std::string& line);

private:
    std::ofstream file;
};