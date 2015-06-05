//
// Created by Kevin Gori on 27/05/15.
//

#include <algorithm>
#include <pll/pll.h>
#include "PLL.h"
std::string PLL::get_tree() {
    pllTreeToNewick(tr->tree_string, tr.get(), partitions,
    tr->start->back, PLL_TRUE, PLL_TRUE,
    PLL_FALSE, PLL_FALSE, PLL_FALSE,
    0, PLL_FALSE, PLL_FALSE);
    std::string tree(tr->tree_string);
    tree.erase(std::remove(tree.begin(), tree.end(), '\n'), tree.end());
    return tree;
}

double PLL::get_likelihood() {
    pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
    return tr->likelihood;
}

int PLL::get_number_of_partitions() {
    return partitions->numberOfPartitions;
}

void PLL::optimise(bool rates, bool freqs, bool alphas, bool branches, double epsilon, bool verbose) {
    if (!rates && !freqs && !alphas && !branches) return;
    int i = 0;
    double loop_start_lnl;
    double loop_end_lnl;
    pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
    for (;;) {
        i++;
        loop_start_lnl = tr->likelihood;
        if (verbose) std::cerr << "  iter " << i << " current lnl = " << loop_start_lnl << std::endl;

        if (rates) {
            pllOptRatesGeneric(tr.get(), partitions, epsilon, partitions->rateList);
            pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
            if (verbose) std::cerr << "    rates:  " << tr->likelihood << std::endl;
        }

        if (branches) {
            pllOptimizeBranchLengths(tr.get(), partitions, 32);
            pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
            if (verbose) std::cerr << "    brlen1: " << tr->likelihood << std::endl;
        }

        if (freqs) {
            pllOptBaseFreqs(tr.get(), partitions, epsilon, partitions->freqList);
            pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
            if (verbose) std::cerr << "    freqs:  " << tr->likelihood << std::endl;
        }

        if (branches) {
            pllOptimizeBranchLengths(tr.get(), partitions, 32);
            pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
            if (verbose) std::cerr << "    brlen2: " << tr->likelihood << std::endl;
        }

        if (alphas) {
            pllOptAlphasGeneric (tr.get(), partitions, epsilon, partitions->alphaList);
            pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
            if (verbose) std::cerr << "    alphas: " << tr->likelihood << std::endl;
        }

        if (branches) {
            pllOptimizeBranchLengths(tr.get(), partitions, 32);
            pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
            if (verbose) std::cerr << "    brlen3: " << tr->likelihood << std::endl;
        }

        loop_end_lnl = tr->likelihood;
        if(loop_end_lnl - loop_start_lnl < 0) {
            std::cerr << loop_end_lnl << " " << loop_start_lnl << std::endl;
            std::cerr << "Difference: " << loop_end_lnl - loop_start_lnl << std::endl;
            break;
        }

        if (loop_end_lnl - loop_start_lnl <= tr->likelihoodEpsilon) {
            if (verbose) {
                std::cerr << "loop_start_lnl = " << loop_start_lnl << std::endl
                          << "loop_end_lnl   = " << loop_end_lnl << std::endl
                          << "END" << std::endl;
            }
            break;
        }
    }
}

const int PLL::get_number_of_partitions() const {
    return partitions->numberOfPartitions;
}

void PLL::tree_search(bool optimise_model) {
    pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
    int pll_bool = optimise_model ? PLL_TRUE : PLL_FALSE;
    pllRaxmlSearchAlgorithm(tr.get(), partitions, pll_bool);
}

void PLL::adjustAlignmentLength(partitionList *partitions, pllAlignmentData *alignment) {
    // Do the adjustment of alignment length in case partition does not cover all sites
    int usedAlLength = 0;
    for (int i = 0; i < partitions->numberOfPartitions; ++i) {
        usedAlLength += partitions->partitionData[i]->width;
    }
    alignment->sequenceLength = usedAlLength;
}

void PLL::set_tree(const std::string& nwk) {
    newickUPtr newick(pllNewickParseString(nwk.c_str()), NewickDeleter());

    if (!newick) {
        throw std::runtime_error("pllNewickParseString returned a null pointer!");
    }

    if (!pllValidateNewick(newick.get())) {
        throw std::runtime_error("Invalid tree!");
    }

    pllTreeInitTopologyNewick(tr.get(), newick.get(), PLL_FALSE);
    pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
}

double PLL::get_alpha(int partition) {
    return pllGetAlpha(partitions, partition);
}

void PLL::set_alpha(double alpha, int partition, bool optimisable) {
    pllSetFixedAlpha(alpha, partition, partitions, tr.get());
    set_optimisable_alpha(partition, optimisable);
}

void PLL::set_optimisable_alpha(int partition, bool optimisable) {
    int pll_bool = optimisable ? PLL_TRUE : PLL_FALSE;
    partitions->partitionData[partition]->optimizeAlphaParameter = pll_bool;
    partitions->dirty = PLL_TRUE;
    pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
}

std::vector<double> PLL::get_frequencies(int partition) {
    std::vector<double> freqs_vec;
    int num_states = partitions->partitionData[partition]->states;
    for (int j = 0; j < num_states; ++j) {
        freqs_vec.push_back(partitions->partitionData[partition]->frequencies[j]);
    }
    return freqs_vec;
}

void PLL::set_frequencies(std::vector<double> freqs, int partition, bool optimisable) {

    double s = 0;
    for (double f : freqs) s += f;
    double diff = 1 - s;
    if (diff < 0) diff = -diff;

    if (diff > EPS) {
        throw std::invalid_argument("Not setting frequencies: Frequencies do not sum to 1");
    }
    size_t num_states = partitions->partitionData[partition]->states;
    if (freqs.size() != num_states) {
        std::ostringstream msg;
        msg << "Frequencies vector is the wrong length. Should be " << num_states;
        throw std::invalid_argument(msg.str());
    }
    set_optimisable_frequencies(partition, true); // frequencies only updated if optimisable flag is true
    pllSetFixedBaseFrequencies(&(freqs[0]), num_states, partition, partitions, tr.get());
    set_optimisable_frequencies(partition, optimisable);
}

void PLL::set_optimisable_frequencies(int partition, bool optimisable) {
    int pll_bool = optimisable ? PLL_TRUE : PLL_FALSE;
    partitions->partitionData[partition]->optimizeBaseFrequencies = pll_bool;
    partitions->dirty = PLL_TRUE;
    pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
}