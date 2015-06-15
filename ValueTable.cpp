//
// Created by Kevin Gori on 15/06/15.
//

#include "ValueTable.h"

ValueTable::ValueTable(unsigned rows, unsigned cols) {
    table = std::vector<std::vector<double>>(rows, std::vector<double>(cols, 0));
    nrow = rows;
    ncol = cols;
}

double ValueTable::get(unsigned r, unsigned c) {
    return table[r][c];
}

std::vector<double> ValueTable::colsum() {
    std::vector<double> s(ncol, 0);
    for (int c = 0; c < ncol; ++c) {
        for (int r = 0; r < nrow; ++r) {
            s[c] += get(r, c);
        }
    }
    return s;
}

std::vector<double> ValueTable::rowsum() {
    std::vector<double> s(nrow, 0);
    for (int r = 0; r < nrow; ++r) {
        for (int c = 0; c < ncol; ++c) {
            s[r] += get(r, c);
        }
    }
    return s;
}

double ValueTable::sum() {
    double s;
    for (int r = 0; r < nrow; ++r) {
        for (int c = 0; c < ncol; ++c) {
            s += get(r, c);
        }
    }
    return s;
}

void ValueTable::set(unsigned r, unsigned c, double val) {
    table[r][c] = val;
}

std::vector<double> ValueTable::rowmean() {
    auto rs = rowsum();
    for (double& val : rs) {
        val /= ncol;
    }
    return rs;
}

std::vector<double> ValueTable::colmean() {
    auto cs = colsum();
    for (double& val : cs) {
        val /= nrow;
    }
    return cs;
}

double ValueTable::mean() {
    double s = sum();
    return s / (nrow * ncol);
}

void ValueTable::print() {
    for (int r=0; r < nrow; ++r) {
        for (int c=0; c < ncol; ++c) {
            std::cout << this->get(r, c) << ' ';
        }
        std::cout << std::endl;
    }
}

std::vector<double> ValueTable::rowmax() {
    std::vector<double> result;
    for (auto& row : table) {
        double max = row[0];
        for (size_t i = 1; i < ncol; ++i) {
            if (row[i] > max) max = row[i];
        }
        result.push_back(max);
    }
    return result;
}
