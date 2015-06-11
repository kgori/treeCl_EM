//
// Created by Kevin Gori on 11/06/15.
//

#ifndef TREECL_EM_OPTIMISER_H
#define TREECL_EM_OPTIMISER_H

#include <vector>

class Optimiser {
public:
    Optimiser();
    void eStep();
    void cStep();
    void mStep();
    std::size_t get_number_of_groups();
    std::vector<int> get_within_group_index(int numgroups);
private:
    std::vector<std::size_t> assignment;
};


#endif //TREECL_EM_OPTIMISER_H
