#include "Matrix.hpp"

Matrix::Matrix(std::vector<std::vector<double>> data) : data(data) {
    rows = data.size();
    cols = data[0].size();
}

Matrix::Matrix(size_t n, size_t m) : rows(n), cols(m) {
    data = std::vector<std::vector<double>>(n, std::vector<double>(m, 0.0));
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
    return (*this) * (1 / x);
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