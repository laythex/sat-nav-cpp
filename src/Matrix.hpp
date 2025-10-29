#pragma once

#include <cmath>
#include <vector>
#include <stdexcept>
#include <iostream>

class Matrix {
public:
    Matrix(std::vector<std::vector<double>> data);
    Matrix(size_t n, size_t m);

    size_t get_rows() const;
    size_t get_cols() const;

    double operator()(size_t i, size_t j) const;
    double& at(size_t i, size_t j);

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(double x) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator/(double x) const;

    Matrix transpose() const;
    Matrix inverse() const;

private:
    size_t rows, cols;
    std::vector<std::vector<double>> data;
};

std::ostream& operator<<(std::ostream& os, const Matrix& matrix);
