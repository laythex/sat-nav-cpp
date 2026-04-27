#include "LinAlg.hpp"

Matrix::Matrix(std::vector<std::vector<double>> data) : data_(data) {
    rows_ = data.size();
    cols_ = data[0].size();
}

Matrix::Matrix(size_t n, size_t m, double x) : rows_(n), cols_(m) {
    data_ = std::vector<std::vector<double>>(n, std::vector<double>(m, x));
}

size_t Matrix::get_rows() const {
    return rows_;
}

size_t Matrix::get_cols() const {
    return cols_;
}

double Matrix::operator()(size_t i, size_t j) const {
    return data_[i][j];
}

double& Matrix::at(size_t i, size_t j) {
    return data_[i][j];
}

std::vector<double> Matrix::operator()(size_t i) const {
    return data_[i];
}

std::vector<double>& Matrix::at(size_t i) {
    return data_[i];
}

void Matrix::insert(const Matrix& m, size_t i, size_t j) {
    for (size_t p = 0; p < m.rows_; ++p) {
        for (size_t q = 0; q < m.cols_; ++q) {
            at(i + p, j + q) = m(p, q);
        }
    }
}

Matrix Matrix::operator+(const Matrix& other) const {
    Matrix res = Matrix(rows_, cols_);

    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < cols_; ++j) {
            res.at(i, j) = operator()(i, j) + other(i, j);
        }
    }

    return res;
}

Matrix Matrix::operator-() const {
    return (*this) * (-1);
}

Matrix Matrix::operator-(const Matrix& other) const {
    return (*this) + other * (-1.0);
}

Matrix Matrix::operator*(double x) const {
    Matrix res = Matrix(rows_, cols_);

    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < cols_; ++j) {
            res.at(i, j) = operator()(i, j) * x;
        }
    }

    return res;
}

std::vector<double> Matrix::operator*(const std::vector<double>& a) const {
    std::vector res(rows_, 0.0);

    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < cols_; ++j) {
            res[i] += operator()(i, j) * a[j];
        }
    }

    return res;
}

Matrix Matrix::operator*(const Matrix& other) const {
    Matrix res = Matrix(rows_, other.cols_);

    for (size_t i = 0; i < rows_; ++i) {
        for (size_t k = 0; k < cols_; ++k) {
            for (size_t j = 0; j < other.cols_; ++j) {
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
    Matrix res = Matrix(cols_, rows_);

    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < cols_; ++j) {
            res.at(j, i) = operator()(i, j);
        }
    }

    return res;
}

Matrix Matrix::inverse() const {
    double singular_tol = 1e-13;

    Matrix aug(rows_, rows_ * 2);

    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < rows_; ++j) {
            aug.at(i, j) = operator()(i, j);
        }
        aug.at(i, rows_ + i) = 1.0;
    }

    double max_pivot = 0.0;
    double min_pivot = std::numeric_limits<double>::infinity();

    for (size_t p = 0; p < rows_; ++p) {

        size_t max_row = p;
        for (size_t i = p + 1; i < rows_; ++i) {
            if (std::abs(aug(i, p)) > std::abs(aug(max_row, p))) {
                max_row = i;
            }
        }

        for (size_t i = 0; i < rows_ * 2; ++i) {
            double tmp = aug(p, i);
            aug.at(p, i) = aug(max_row, i);
            aug.at(max_row, i) = tmp;
        }

        double pivot = std::abs(aug(p, p));
        if (pivot > max_pivot) max_pivot = pivot;
        if (pivot < min_pivot) min_pivot = pivot;
        if (pivot <= singular_tol * max_pivot) {
            throw std::runtime_error("Singular matrix");
        }

        for (size_t i = p + 1; i < rows_; ++i) {
            double factor = aug(i, p) / aug(p, p);
            for (size_t j = p; j < rows_ * 2; ++j) {
                aug.at(i, j) -= factor * aug(p, j);
            }
        }
    }

    double rough_cond = max_pivot / min_pivot;
    if (rough_cond > 1e8) {
        throw std::runtime_error("Ill-conditioned matrix");
    }

    for (size_t p1 = rows_; p1 > 0; p1--) {
        size_t p = p1 - 1;

        double pivot = aug(p, p);
        for (size_t j = 0; j < rows_ * 2; ++j) {
            aug.at(p, j) /= pivot;
        }

        for (size_t i = 0; i < p; ++i) {
            double factor = aug(i, p);
            for (size_t j = 0; j < rows_ * 2; ++j) {
                aug.at(i, j) -= factor * aug(p, j);
            }
        }
    }

    Matrix res(rows_, rows_);

    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < rows_; ++j) {
            res.at(i, j) = aug(i, rows_ + j);
        }
    }

    return res;
}

double Matrix::trace() const {
    double res = 0;

    for (size_t i = 0; i < rows_; ++i) {
        res += operator()(i, i);
    }

    return res;
}

double Matrix::norm() const {
    double res = 0;

    for (size_t i = 0; i < rows_; ++i) {
        double row_abs_sum = 0;
        for (size_t j = 0; j < cols_; ++j) {
            row_abs_sum += std::abs(operator()(i, j));
        }
        if (row_abs_sum > res) {
            res = row_abs_sum;
        }
    }

    return res;
}

std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    for (size_t i = 0; i < matrix.get_rows(); ++i) {
        for (size_t j = 0; j < matrix.get_cols(); ++j) {
            double x = std::abs(matrix(i, j)) > 1e-13 ? matrix(i, j) : 0;
            os << x << " ";
        }
        os << '\n';
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<double>& a) {
    for (size_t i = 0; i < a.size(); ++i) {
        os << a[i] << " ";
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<unsigned>& a) {
    for (size_t i = 0; i < a.size(); ++i) {
        double x = a[i];
        os << x << " ";
    }

    return os;
}

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> res(a.size());

    for (size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] + b[i];
    }

    return res;
}

std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b) {
    return a + b * (-1.0);
}

std::vector<double> operator*(const std::vector<double>& a, double x) {
    std::vector<double> res(a.size());

    for (size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] * x;
    }

    return res;
}

std::vector<double> operator*(double x, const std::vector<double>& a) {
    return a * x;
}

double operator*(const std::vector<double>& a, const std::vector<double>& b) {
    return dot(a, b);
}

std::vector<double> operator/(const std::vector<double>& a, double x) {
    return a * (1.0 / x);
}

Matrix operator*(double x, const Matrix& matrix) {
    return matrix * x;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double res = 0;

    for (size_t i = 0; i < a.size(); ++i) {
        res += a[i] * b[i];
    }

    return res;
}

Matrix tensor(const std::vector<double>& a, const std::vector<double>& b) {
    Matrix res(a.size(), b.size());

    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            res.at(i, j) = a[i] * b[j];
        }
    }

    return res;
}

double abs(const std::vector<double>& a) {
    return sqrt(dot(a, a));
}

double abs(const Matrix& matrix) {
    return matrix.norm();
}

double angle_between(const std::vector<double>& a, const std::vector<double>& b) {
    return acos(dot(a, b) / abs(a) / abs(b));
}

std::vector<double> normalize(const std::vector<double>& a) {
    return a / abs(a);
}

Matrix zero(size_t r, size_t c) {
    return Matrix(r, c);
}

Matrix identity(size_t r) {
    Matrix res(r, r);

    for (size_t i = 0; i < r; ++i) {
        res.at(i, i) = 1.0;
    }

    return res;
}

Matrix rotation(double angle, Axis axis) {
    Matrix res = identity(3);

    double c = cos(angle);
    double s = sin(angle);

    switch(axis) {
        case Axis::X:
            res.at(1, 1) = c;
            res.at(2, 2) = c;
            res.at(1, 2) = -s;
            res.at(2, 1) = s;
            break;
        case Axis::Y:
            res.at(0, 0) = c;
            res.at(2, 2) = c;
            res.at(0, 2) = s;
            res.at(2, 0) = -s;
            break;
        case Axis::Z:
            res.at(0, 0) = c;
            res.at(1, 1) = c;
            res.at(0, 1) = -s;
            res.at(1, 0) = s;
            break;
    }

    return res;
}

double dist_to_bounds(double x, double lo, double hi) {
    return std::min(x - lo, hi - x);
}

size_t find_max_abs(const std::vector<double>& a) {
    size_t max_at = 0;
    double max_el = 0;

    for (size_t i = 0; i < a.size(); ++i) {
        double el = std::abs(a[i]);
        if (el > max_el) {
            max_el = el;
            max_at = i;
        }
    }

    return max_at;
}