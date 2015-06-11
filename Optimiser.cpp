//
// Created by Kevin Gori on 11/06/15.
//

#include "Optimiser.h"
#include <algorithm>

size_t Optimiser::get_number_of_groups() {
    auto max_elem = std::max_element(assignment.begin(), assignment.end());
    return 1 + assignment[std::distance(assignment.begin(), max_elem)];
}

std::vector<int> Optimiser::get_within_group_index(int numgroups) {
    std::vector<int> grpix(numgroups, 0);
    std::vector<int> result(assignment.size(), 0);
    for (int i = 0; i < assignment.size(); ++i) {
        int grp = assignment[i];
        result[i] = grpix[grp]++;
    }
    return result;
}
