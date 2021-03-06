//
// Created by Kevin Gori on 15/06/15.
//

#ifndef TREECL_EM_UTILS_H
#define TREECL_EM_UTILS_H

#include <sstream>
#include <vector>
#include "memory_management.h"

namespace utils {
    bool is_file(std::string filename);
    std::vector<std::string> readlines(const std::string& filename);
    alignmentUPtr parse_alignment_file(std::string path);
    queueUPtr parse_partitions(std::string partitions);
    double logsumexp(const std::vector<double>& nums);
    std::vector<double> scale_by_sum(const std::vector<double>& nums);
    size_t random_select(const std::vector<double>& probs);

    template<typename T>
    std::vector<T> csum(const std::vector<T>& row);

    template <typename Iter>
    std::string stringjoin(Iter it, Iter end, char delim = '\n') {
        std::ostringstream oss;
        if (it == end) return oss.str();
        --end;
        for (;it != end; ++it) {
            oss << *it << delim;
        }
        oss << *(it);
        return oss.str();
    }

    template<typename Iter>
    void print_container(Iter it, Iter end, char delim = ' ') {
        if (it == end) return;
        --end;
        for (;it != end; ++it) {
            std::cout << *it << " ";
        }
        std::cout << *(it) << std::endl;
    }
}

#endif //TREECL_EM_UTILS_H
