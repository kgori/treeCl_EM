//
// Created by Kevin Gori on 27/05/15.
//

#ifndef TREECL_EM_PLL_H
#define TREECL_EM_PLL_H

#include "memory_management.h"

class PLL{

public:
    PLL(pllInstanceAttr& attr, pllQueue* queue, pllAlignmentData* alignment) {
        tr = instanceUPtr(pllCreateInstance(&attr));
        partitions = pllPartitionsCommit(queue, alignment);
        pllAlignmentRemoveDups(alignment, partitions);
        pllTreeInitTopologyRandom(tr.get(), alignment->sequenceCount, alignment->sequenceLabels);
        pllLoadAlignment(tr.get(), alignment, partitions);
        pllComputeRandomizedStepwiseAdditionParsimonyTree(tr.get(), partitions);
        pllInitModel(tr.get(), partitions);
    }

    PLL(PLL&& rhs) = delete;

    PLL& operator=(PLL&& rhs) = delete;

    ~PLL() {
        if (partitions) {
            pllPartitionsDestroy(tr.get(), &partitions);
            std::cout << "Destroyed partitions" << std::endl;
        }
    }

    void optimise(bool rates, bool freqs, bool alphas, bool branches, double epsilon=0.0001);

    partitionList* partitions;
    instanceUPtr tr;
};


#endif //TREECL_EM_PLL_H
