//
// Created by Kevin Gori on 27/05/15.
//

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
