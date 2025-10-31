#include "LinearAlgebra.hpp"

Matrix::Matrix(std::vector<std::vector<double>> data) : data(data) {
    rows = data.size();
    cols = data[0].size();
}

Matrix::Matrix(size_t n, size_t m, double x) : rows(n), cols(m) {
    data = std::vector<std::vector<double>>(n, std::vector<double>(m, x));
}

size_t Matrix::get_rows() const {
    return rows;
}

size_t Matrix::get_cols() const {
    return cols;
}

double Matrix::operator()(size_t i, size_t j) const {
    return data[i][j];
}

double& Matrix::at(size_t i, size_t j) {
    return data[i][j];
}

Matrix Matrix::operator+(const Matrix& other) const {
    Matrix res = Matrix(rows, cols);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            res.at(i, j) = operator()(i, j) + other(i, j);
        }
    }

    return res;
}

Matrix Matrix::operator-(const Matrix& other) const {
    return (*this) + other * (-1.0);
}

Matrix Matrix::operator*(double x) const {
    Matrix res = Matrix(rows, cols);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            res.at(i, j) = operator()(i, j) * x;
        }
    }

    return res;
}

std::vector<double> Matrix::operator*(const std::vector<double>& a) const {
    std::vector res(rows, 0.0);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            res[i] += operator()(i, j) * a[j];
        }
    }

    return res;
}

Matrix Matrix::operator*(const Matrix& other) const {
    Matrix res = Matrix(rows, other.cols);

    for (size_t i = 0; i < rows; i++) {
        for (size_t k = 0; k < cols; k++) {
            for (size_t j = 0; j < other.cols; j++) {
                res.at(i, j) += operator()(i, k) * other(k, j);
            }
        }
    }
    
    return res;
}

Matrix Matrix::operator/(double x) const {
    return (*this) * (1.0 / x);
}

Matrix Matrix::transpose() const {
    Matrix res = Matrix(cols, rows);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            res.at(j, i) = operator()(i, j);
        }
    }

    return res;
}

Matrix Matrix::inverse() const {
    Matrix aug(rows, rows * 2);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < rows; j++) {
            aug.at(i, j) = operator()(i, j);
        }
        aug.at(i, rows + i) = 1.0;
    }

    for (size_t p = 0; p < rows; p++) {

        size_t max_row = p;
        for (size_t i = p + 1; i < rows; i++) {
            if (std::abs(aug(i, p)) > std::abs(aug(max_row, p))) {
                max_row = i;
            }
        }

        for (size_t i = 0; i < rows * 2; i++) {
            double tmp = aug(p, i);
            aug.at(p, i) = aug(max_row, i);
            aug.at(max_row, i) = tmp;
        }

        if (std::abs(aug(p, p)) < 1e-13) {
            throw std::runtime_error("Singular matrix");
        }

        for (size_t i = p + 1; i < rows; i++) {
            double factor = aug(i, p) / aug(p, p);
            for (size_t j = p; j < rows * 2; j++) {
                aug.at(i, j) -= factor * aug(p, j);
            }
        }
    }

    for (size_t p1 = rows; p1 > 0; p1--) {
        size_t p = p1 - 1;

        double pivot = aug(p, p);
        for (size_t j = 0; j < rows * 2; j++) {
            aug.at(p, j) /= pivot;
        }

        for (size_t i = 0; i < p; i++) {
            double factor = aug(i, p);
            for (size_t j = 0; j < rows * 2; j++) {
                aug.at(i, j) -= factor * aug(p, j);
            }
        }
    }

    Matrix res(rows, rows);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < rows; j++) {
            res.at(i, j) = aug(i, rows + j);
        }
    }

    return res;
}

std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    for (size_t i = 0; i < matrix.get_rows(); i++) {
        for (size_t j = 0; j < matrix.get_cols(); j++) {
            double x = std::abs(matrix(i, j)) > 1e-13 ? matrix(i, j) : 0;
            os << x << " ";
        }
        os << '\n';
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<double>& a) {
    for (size_t i = 0; i < a.size(); i++) {
        double x = abs(a[i]) > 1e-13 ? a[i] : 0;
        os << x << " ";
    }

    return os;
}


std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> res(a.size());

    for (size_t i = 0; i < a.size(); i++) {
        res[i] = a[i] + b[i];
    }

    return res;
}

std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b) {
    return a + b * (-1.0);
}

std::vector<double> operator*(const std::vector<double>& a, double x) {
    std::vector<double> res(a.size());

    for (size_t i = 0; i < a.size(); i++) {
        res[i] = a[i] * x;
    }

    return res;
}

std::vector<double> operator/(const std::vector<double>& a, double x) {
    return a * (1.0 / x);
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double res = 0;

    for (size_t i = 0; i < a.size(); i++) {
        res += a[i] * b[i];
    }

    return res;
}

double abs(const std::vector<double>& a) {
    return sqrt(dot(a, a));
}

Matrix identity(size_t s) {
    Matrix res(s, s);

    for (size_t i = 0; i < s; i++) {
        res.at(i, i) = 1.0;
    }

    return res;
}

Matrix rotation(double angle, char axis) {
    Matrix res = identity(3);

    double c = cos(angle);
    double s = sin(angle);

    switch(axis) {
        case 'x':
            res.at(1, 1) = c;
            res.at(2, 2) = c;
            res.at(1, 2) = -s;
            res.at(2, 1) = s;
            break;
        case 'y':
            res.at(0, 0) = c;
            res.at(2, 2) = c;
            res.at(0, 2) = s;
            res.at(2, 0) = -s;
            break;
        case 'z':
            res.at(0, 0) = c;
            res.at(1, 1) = c;
            res.at(0, 1) = -s;
            res.at(1, 0) = s;
            break;
    }

    return res;
}