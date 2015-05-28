//
// Created by Kevin Gori on 27/05/15.
//

#ifndef TREECL_EM_PLL_H
#define TREECL_EM_PLL_H

#include <pll/pll.h>
#include "memory_management.h"
typedef pInfo *pInfoPtr;
class PLL{

public:
    PLL(pllInstanceAttr& attr, pllQueue* queue, pllAlignmentData* alignment) {
        tr = instanceUPtr(pllCreateInstance(&attr));
        partitions = pllPartitionsCommit(queue, alignment);
        pllAlignmentRemoveDups(alignment, partitions);

        // Here do the adjustment of alignment length in case partition does not cover all sites
        adjustAlignmentLength(partitions, alignment);

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

        // Here do the adjustment of alignment length in case partition does not cover all sites
        adjustAlignmentLength(partitions, alignment);

        pllTreeInitTopologyNewick(tr.get(), newick, PLL_FALSE);
        if (!pllLoadAlignment(tr.get(), alignment, partitions)) {
            std::cout << "Problem" << std::endl;
        }
        pllComputeRandomizedStepwiseAdditionParsimonyTree(tr.get(), partitions);
        pllInitModel(tr.get(), partitions);
    }

    // Move constructor
    PLL (PLL&& other) noexcept : partitions(std::move(other.partitions)), tr(std::move(other.tr)) {
        std::cout << "Move constructed PLL" << std::endl;
    }

    /** Move assignment operator */
    PLL& operator= (PLL&& other) noexcept {
        // simplified move-constructor that also protects against move-to-self.
        partitions = std::move(other.partitions);
        tr = std::move(other.tr);
        std::cout << "Move assigned PLL" << std::endl;
        return *this;
    }

    // Delete copy assignment
    PLL& operator=(const PLL& other) = delete;

    // Delete copy constructor
    PLL(const PLL& other) = delete;

    // Destructor
    ~PLL() noexcept {
        if (partitions) {
            pllPartitionsDestroy(tr.get(), &partitions);
            std::cout << "Destroyed partitions" << std::endl;
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

    void tree_search(bool optimise_model);

public:
    partitionList* partitions;
    instanceUPtr tr;
    void adjustAlignmentLength(partitionList* partitions, pllAlignmentData* alignment);
};

typedef std::unique_ptr<PLL> PLLUPtr;
typedef PLL *PLLPtr;

#endif //TREECL_EM_PLL_H
