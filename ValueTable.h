//
// Created by Kevin Gori on 15/06/15.
//

#ifndef TREECL_EM_VALTABLE_H
#define TREECL_EM_VALTABLE_H

#include <iostream>
#include <vector>

class ValueTable {
    std::vector<std::vector<double>> table;
    unsigned nrow;
    unsigned ncol;

public:
    ValueTable(unsigned rows, unsigned cols);
    unsigned nrows() { return nrow; }
    unsigned ncols() { return ncol; }
    double get(unsigned r, unsigned c);
    void set(unsigned r, unsigned c, double val);
    std::vector<double> rowsum();
    std::vector<double> colsum();
    double sum();
    std::vector<double> rowmean();
    std::vector<double> colmean();
    double mean();
    std::vector<double> rowmax();
    const std::vector<std::vector<double>>& get_table() { return table; };
    void print();
};


#endif //TREECL_EM_VALTABLE_H
