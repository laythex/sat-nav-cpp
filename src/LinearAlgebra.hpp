#pragma once

#include <cmath>
#include <vector>
#include <stdexcept>
#include <iostream>

class Matrix {
public:
    Matrix(std::vector<std::vector<double>> data);
    Matrix(size_t n, size_t m, double x = 0.0);

    size_t get_rows() const;
    size_t get_cols() const;

    double operator()(size_t i, size_t j) const;
    double& at(size_t i, size_t j);

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(double x) const;
    std::vector<double> operator*(const std::vector<double>& a) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator/(double x) const;

    Matrix transpose() const;
    Matrix inverse() const;
    double trace() const;

private:
    size_t rows, cols;
    std::vector<std::vector<double>> data;
};

std::ostream& operator<<(std::ostream& os, const Matrix& matrix);
std::ostream& operator<<(std::ostream& os, const std::vector<double>& a);

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> operator*(const std::vector<double>& a, double x);
std::vector<double> operator/(const std::vector<double>& a, double x);

double dot(const std::vector<double>& a, const std::vector<double>& b);
double abs(const std::vector<double>& a);

Matrix identity(size_t s);
Matrix rotation(double angle, char axis);