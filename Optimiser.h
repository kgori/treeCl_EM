//
// Created by Kevin Gori on 11/06/15.
//

#ifndef TREECL_EM_OPTIMISER_H
#define TREECL_EM_OPTIMISER_H

#include <map>
#include <string>
#include <vector>
#include <limits>
#include "memory_management.h"
#include "PLL.h"
#include "utils.h"

const double UNLIKELY = std::numeric_limits<double>::lowest();

// Parameters that belong to the locus
struct perlocus {
    double alpha;
    std::vector<double> freqs;
    std::vector<double> rates;
    double likelihood;
};

//Parameters that belong to the group
struct pergroup {
    std::string tree;
    double likelihood;
};

//PLL result
struct pllresult {
    std::vector<perlocus> locus_params;
    std::string tree;
    double likelihood;
};

enum class Schedule {
    NO_SEARCH,      // Only optimise branch lengths
    PARAM_SEARCH,   // Optimise branch lengths and model parameters, but keep input tree
    TREE_SEARCH,    // Optimise branch lengths and do tree search, don't optimise model parameters
    FULL_SEARCH     // Optimise everything and do tree search
}; // Optimisation schedule

enum class Classifier {
    MAP,    // Assign to maximum aposteriori probability group
    IMPUTE, // Stochastic assignment due to posterior probability distribution
    ANNEAL, // Simulated annealing
}; // Classification criterion



class Optimiser {
public:
    Optimiser(const std::string alignment, const std::vector<std::string>& partitions, attrSPtr attr) :
        alignment(alignment), partitions(partitions), attr(attr), nLoci(partitions.size()) {
        parameters.resize(nLoci);
    };
    void set_assignment(const std::vector<int>& a);
    void set_assignment(int nGroups);
    void set_qs(const std::vector<int>& a);
    void eStep();
    void cStep();
    void mStep();
    pllresult doOpt(PLLUPtr&& pll, Schedule schedule);
    //void doIteration();
    double get_likelihood() { return likelihood; };
    const std::vector<int>& get_assignment() { return assignment; };
    int get_number_of_groups(const std::vector<int>& a);
    void index(const std::vector<int>& a);
    std::vector<int> make_random_assignment();
    std::vector<double> get_proportions(int pseudocount=1);
    pllresult get_parameters(PLLUPtr&& pll);

private:
    unsigned nGroups;
    unsigned nLoci;
    const std::string alignment;
    std::vector<int> assignment;
    std::map<int, std::vector<int>> indexmap;
    std::vector<std::unique_ptr<std::stringstream>> qs;
    const std::vector<std::string>& partitions;
    double likelihood = UNLIKELY;
    queueUPtr q;
    alignmentUPtr al;
    attrSPtr attr;
    Schedule schedule = Schedule::NO_SEARCH;
    Classifier classifier = Classifier::MAP;
    std::vector<perlocus> parameters;
    std::vector<pergroup> trees;
    std::vector<double> proportions;
    bool have_parameters = false;
};


#endif //TREECL_EM_OPTIMISER_H
