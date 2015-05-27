//
// Created by Kevin Gori on 27/05/15.
//

#ifndef TREECL_EM_PLL_H
#define TREECL_EM_PLL_H

#include "memory_management.h"
typedef pInfo *pInfoPtr;
class PLL{

public:
    PLL(pllInstanceAttr& attr, pllQueue* queue, pllAlignmentData* alignment) {
        tr = instanceUPtr(pllCreateInstance(&attr));
        partitions = pllPartitionsCommit(queue, alignment);
        pllAlignmentRemoveDups(alignment, partitions);
        pllTreeInitTopologyRandom(tr.get(), alignment->sequenceCount, alignment->sequenceLabels);
        if (!pllLoadAlignment(tr.get(), alignment, partitions)) {
            std::cout << "Problem" << std::endl;
        }
        pllComputeRandomizedStepwiseAdditionParsimonyTree(tr.get(), partitions);
        pllInitModel(tr.get(), partitions);
        pllTreeToNewick(tr->tree_string, tr.get(), partitions,
                        tr->start->back, PLL_TRUE, PLL_TRUE,
                        PLL_FALSE, PLL_FALSE, PLL_FALSE,
                        0, PLL_FALSE, PLL_FALSE);
        newickUPtr newick = newickUPtr(pllNewickParseString(tr->tree_string));
        pllTreeInitTopologyNewick(tr.get(), newick.get(), PLL_TRUE);
    }

    PLL(pllInstanceAttr& attr, pllQueue* queue, pllNewickTree* newick, pllAlignmentData* alignment) {
        tr = instanceUPtr(pllCreateInstance(&attr));
        partitions = pllPartitionsCommit(queue, alignment);
        pllAlignmentRemoveDups(alignment, partitions);
        pllTreeInitTopologyNewick(tr.get(), newick, PLL_FALSE);
        if (!pllLoadAlignment(tr.get(), alignment, partitions)) {
            std::cout << "Problem" << std::endl;
        }
        pllComputeRandomizedStepwiseAdditionParsimonyTree(tr.get(), partitions);
        pllInitModel(tr.get(), partitions);
    }

    PLL(PLL&& rhs) = delete;

    PLL& operator=(PLL&& rhs) = delete;

    ~PLL() {
        if (partitions) {
            pllPartitionsDestroy(tr.get(), &partitions);
            //std::cout << "Destroyed partitions" << std::endl;
        }
    }

    void optimise(bool rates, bool freqs, bool alphas, bool branches, double epsilon=0.0001, bool verbose=false);
    std::string get_tree();
    double get_likelihood();
    int get_number_of_partitions();
    const int get_number_of_partitions() const;
    pInfoPtr& operator[](int i)  {
        if (i > get_number_of_partitions()) throw std::exception();
        return partitions->partitionData[i];
    };
    const pInfoPtr& operator[](int i) const  {
        if (i > get_number_of_partitions()) throw std::exception();
        return partitions->partitionData[i];
    };

private:
    partitionList* partitions;
    instanceUPtr tr;
};


#endif //TREECL_EM_PLL_H
