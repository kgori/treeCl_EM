//
// Created by Kevin Gori on 27/05/15.
//

#ifndef TREECL_EM_PLL_H
#define TREECL_EM_PLL_H

extern "C" {
    #include <pll/pll.h>
}
#include <sstream>
#include <stdexcept>
#include "memory_management.h"
#include "threadpool.h"

#define EPS 1e-6

using pInfoPtr = pInfo*;
//typedef pInfo *pInfoPtr;
class PLL{

public:
    PLL(pllInstanceAttr& attr, pllQueue* queue, pllAlignmentData* alignment) {
        tr = instanceUPtr(pllCreateInstance(&attr), InstanceDeleter());
        partitions = pllPartitionsCommit(queue, alignment);
        pllAlignmentRemoveDups(alignment, partitions);

        // Here do the adjustment of alignment length in case partition does not cover all sites
        // THIS IS A HACK TO TRY TO MAKE THREADING WORK WITH INCOMPLETE PARTITION COVERAGE
        // Given the memory issues with multiple PLL instances all using multithreading,
        // this is the least of my worries
        adjustAlignmentLength(partitions, alignment);

        pllTreeInitTopologyRandom(tr.get(), alignment->sequenceCount, alignment->sequenceLabels);
        if (!pllLoadAlignment(tr.get(), alignment, partitions)) {
            std::cout << "Problem" << std::endl;
        }
        pllComputeRandomizedStepwiseAdditionParsimonyTree(tr.get(), partitions);
        pllInitModel(tr.get(), partitions);

        // Here reinitialise the tree because, for some reason, pllOptimizeBranchLengths doesn't
        // work with random trees or parsimony trees. Maybe it's the tiny branches, who knows?
        pllTreeToNewick(tr->tree_string, tr.get(), partitions,
                        tr->start->back, PLL_TRUE, PLL_TRUE,
                        PLL_FALSE, PLL_FALSE, PLL_FALSE,
                        0, PLL_FALSE, PLL_FALSE);
        newickUPtr newick = newickUPtr(pllNewickParseString(tr->tree_string));
        pllTreeInitTopologyNewick(tr.get(), newick.get(), PLL_TRUE);
    }

    PLL(pllInstanceAttr& attr, pllQueue* queue, pllNewickTree* newick, pllAlignmentData* alignment) {
        tr = instanceUPtr(pllCreateInstance(&attr), InstanceDeleter());
        partitions = pllPartitionsCommit(queue, alignment);
        pllAlignmentRemoveDups(alignment, partitions);

        // Here do the adjustment of alignment length in case partition does not cover all sites
        adjustAlignmentLength(partitions, alignment);

        pllTreeInitTopologyNewick(tr.get(), newick, PLL_FALSE);
        if (!pllLoadAlignment(tr.get(), alignment, partitions)) {
            throw std::runtime_error("Error loading alignment");
        }
        pllInitModel(tr.get(), partitions);
    }

    // Move constructor
    PLL (PLL&& other) noexcept : partitions(other.partitions), tr(std::move(other.tr)) {
        other.partitions = nullptr;
    }

    /** Move assignment operator */
    PLL& operator= (PLL&& other) noexcept {
        // simplified move-constructor that also protects against move-to-self.
        std::swap(partitions, other.partitions);
        std::swap(tr, other.tr);
        return *this;
    }

    // Destructor
    ~PLL() noexcept {
        if (partitions && tr) {
            pllPartitionsDestroy(tr.get(), &partitions);
        }
    }

    // Copy constructor
    PLL(const PLL& other) = delete;

    // Copy assignment
    PLL& operator=(const PLL& other) = delete;

    void optimise(bool rates, bool freqs, bool alphas, bool branches, double epsilon=0.0001, bool verbose=false);
    std::string get_tree();
    double get_likelihood();
    int get_number_of_partitions();
    const int get_number_of_partitions() const;
    pInfoPtr& operator[](int i)  {
        if (i > get_number_of_partitions()) throw std::out_of_range("Index out of bounds");
        return partitions->partitionData[i];
    };
    const pInfoPtr& operator[](int i) const  {
        if (i > get_number_of_partitions()) throw std::out_of_range("Index out of bounds");
        return partitions->partitionData[i];
    };

    unsigned sites();

    void tree_search(bool optimise_model);
    void set_tree(const std::string& nwk);

    double get_alpha(int partition);
    void set_alpha(double alpha, int partition, bool optimisable);
    void set_optimisable_alpha(int partition, bool optimisable);

    std::vector<double> get_frequencies(int partition);
    void set_frequencies(std::vector<double> freqs, int partition, bool optimisable);
    void set_optimisable_frequencies(int partition, bool optimisable);

    std::vector<double> get_rates(int partition);
    void set_rates(std::vector<double> rates, int partition, bool optimisable);
    void set_optimisable_rates(int partition, bool optimisable);

public:
    partitionList* partitions = nullptr;
    instanceUPtr tr;
    void adjustAlignmentLength(partitionList* partitions, pllAlignmentData* alignment);
};

typedef std::unique_ptr<PLL> PLLUPtr;

#endif //TREECL_EM_PLL_H
