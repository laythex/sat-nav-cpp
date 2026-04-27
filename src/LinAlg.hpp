#pragma once

#include <cmath>
#include <vector>
#include <stdexcept>
#include <iostream>

enum class Axis {X, Y, Z};

class Matrix {
public:
    Matrix(std::vector<std::vector<double>> data);
    Matrix(size_t n, size_t m, double x = 0.0);

    size_t get_rows() const;
    size_t get_cols() const;

    double operator()(size_t i, size_t j) const;
    double& at(size_t i, size_t j);

    std::vector<double> operator()(size_t i) const;
    std::vector<double>& at(size_t i);
    void insert(const Matrix& m, size_t i, size_t j);

    Matrix operator+(const Matrix& other) const;
    Matrix operator-() const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(double x) const;
    std::vector<double> operator*(const std::vector<double>& a) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator/(double x) const;

    Matrix transpose() const;
    Matrix inverse() const;
    double trace() const;
    double norm() const;

private:
    size_t rows_, cols_;
    std::vector<std::vector<double>> data_;
};

std::ostream& operator<<(std::ostream& os, const Matrix& matrix);
std::ostream& operator<<(std::ostream& os, const std::vector<double>& a);
std::ostream& operator<<(std::ostream& os, const std::vector<unsigned>& a);

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> operator*(const std::vector<double>& a, double x);
std::vector<double> operator*(double x, const std::vector<double>& a);
double operator*(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> operator/(const std::vector<double>& a, double x);

Matrix operator*(double x, const Matrix& matrix);

double dot(const std::vector<double>& a, const std::vector<double>& b);
Matrix tensor(const std::vector<double>& a, const std::vector<double>& b);

double abs(const std::vector<double>& a);
double abs(const Matrix& matrix);
double angle_between(const std::vector<double>& a, const std::vector<double>& b);

std::vector<double> normalize(const std::vector<double>& a);

Matrix zero(size_t r, size_t c);
Matrix identity(size_t r);
Matrix rotation(double angle, Axis axis);

double dist_to_bounds(double x, double lo, double hi);

size_t find_max_abs(const std::vector<double>& a);

template <typename T>
std::vector<T> mask(const std::vector<T>& a, const std::vector<bool> m);

template <typename T>
std::vector<T> mask(const std::vector<T>& a, const std::vector<bool> m) {
    std::vector<T> res;

    for (size_t i = 0; i < a.size(); i++) {
        if (m[i]) res.push_back(a[i]);
    }

    return res;
}
