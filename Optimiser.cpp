//
// Created by Kevin Gori on 11/06/15.
//

#include "Optimiser.h"
#include "utils.h"
#include "ValueTable.h"
#include <algorithm>
#include <cmath>
#include <random>

int Optimiser::get_number_of_groups(const std::vector<int>& a) {
    auto max_elem = std::max_element(a.begin(), a.end());
    return 1 + *max_elem;
}

void Optimiser::index(const std::vector<int>& a) {
    indexmap.clear();
    for (int i = 0; i < a.size(); ++i) {
        int grp = a[i];
        indexmap[grp].push_back(i);
    }
}

std::vector<int> Optimiser::make_random_assignment() {
    std::vector<int> v;
    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<int> dist(0, nGroups-1);

    for (int i=0; i < nGroups; ++i) {
        v.push_back(i);
    }

    for (int i=nGroups; i < nLoci; ++i ) {
        v.push_back(dist(engine));
    }

    std::shuffle(v.begin(), v.end(), engine);
    return v;
}

std::vector<double> Optimiser::get_proportions(int pseudocount) {
    std::vector<double> props;

    for (int i=0; i<nGroups; ++i) {
        props.push_back(std::count(assignment.begin(), assignment.end(), i));
    }

    for (double& n : props) {
        n = (n + pseudocount) / (nLoci + (pseudocount * nGroups));
    }

    return props;
}

void Optimiser::make_probability_table() {
    std::vector<double> logprobsum;
    for (const auto& row : vtab->get_table()) {
        double lps = utils::logsumexp(row);
        std::cout << "LPS = " << lps << std::endl;
        logprobsum.push_back(lps);
    }

    for (int i=0; i < vtab->nrows(); ++i) {
        for (int j=0; j < vtab->ncols(); ++j) {
            vtab->set(i, j, exp(vtab->get(i, j) - logprobsum[i]));
        }
    }
}

void Optimiser::eStep() {
    for (int i=0; i < nLoci; ++i) {
        // TODO: nLoci reads of the alignment file every eStep? Too many!
        al = utils::parse_alignment_file(alignment);
        q = utils::parse_partitions(partitions[i].c_str());
        PLLUPtr pll = std::make_unique<PLL>(*attr, q.get(), al.get());

        // Load current parameter estimates
        pll->set_alpha(parameters[i].alpha, 0, false);
        pll->set_frequencies(parameters[i].freqs, 0, false);
        pll->set_rates(parameters[i].rates, 0, false);
        for (int j=0; j < nGroups; ++j) {
            auto x=trees.size();
            auto tree = trees[j].tree;
            pll->set_tree(trees[j].tree);
            vtab->set(i, j, pll->get_likelihood() + log(proportions[j]));
        }
    }
    std::cout << "INTERMEDIATE VTAB" << std::endl;
    vtab->print();
    make_probability_table();
}

void Optimiser::cStep() {
    switch (classifier) {
        case Classifier::MAP:
            break;
        case Classifier::IMPUTE:
            break;
        case Classifier::ANNEAL:
            break;
    }
}

void Optimiser::mStep() {
    int index;
    double lnl = 0;

    // Update phylogenetic parameters
    for (int g = 0; g < nGroups; ++g) {
        // TODO: Don't recalculate if a group hasn't changed
        al = utils::parse_alignment_file(alignment);
        std::string qstring = qs[g]->str();
        q = utils::parse_partitions(qstring.c_str());
        PLLUPtr pll = std::make_unique<PLL>(*attr, q.get(), al.get());

        // Load current parameter estimates, if we have any yet (only don't on first iteration)
        if (have_parameters) {
            bool opt = (schedule == Schedule::PARAM_SEARCH || schedule == Schedule::FULL_SEARCH);
            for (int wgi = 0;
                 wgi < indexmap[g].size(); ++wgi) {  // wgi = within group index; wdi = within dataset index
                int wdi = indexmap[g][wgi];
                pll->set_alpha(parameters[wdi].alpha, wgi, opt);
                pll->set_frequencies(parameters[wdi].freqs, wgi, opt);
                pll->set_rates(parameters[wdi].rates, wgi, opt);
            }
            pll->set_tree(trees[g].tree);
        }

        // Optimise
        auto result = doOpt(std::move(pll), schedule);

        // Save parameters
        for (int i=0; i < result.locus_params.size(); ++i) {
            index = indexmap[g][i];
            parameters[index] = result.locus_params[i];
        }
        trees[g].tree = result.tree;
        trees[g].likelihood = result.likelihood;
    };
    have_parameters = true;

    // Update assignment probabilities
    proportions = get_proportions();

    for (auto& tree: trees) {
        lnl += tree.likelihood;
    }
    likelihood = lnl;
}

pllresult Optimiser::doOpt(PLLUPtr&& pll, Schedule schedule) {
    switch(schedule) {
        case Schedule::NO_SEARCH:
            pll->optimise(false, false, false, true, EPS, false);
            break;
        case Schedule::PARAM_SEARCH:
            pll->optimise(true, true, true, true, EPS, false);
            break;
        case Schedule::TREE_SEARCH:
            pll->tree_search(false);
            break;
        case Schedule::FULL_SEARCH:
            pll->tree_search(true);
            break;
    }
    return get_parameters(std::move(pll));
}


void Optimiser::set_assignment(const std::vector<int>& a) {
    nGroups = get_number_of_groups(a);
    trees.clear();
    trees.resize(nGroups);
    assignment = a;
    index(a);
    set_qs(a);
    vtab = std::make_unique<ValueTable>(nLoci, nGroups);
}

void Optimiser::set_assignment(int nGroups) {
    this->nGroups = nGroups;
    auto a = make_random_assignment();
    set_assignment(a);
}

void Optimiser::set_qs(const std::vector<int>& a) {
    qs.clear();
    for (int i = 0; i < nGroups; ++i) {
        qs.push_back(std::make_unique<std::stringstream>());
    }
    for (int j = 0; j < a.size(); ++j) {
        int group_index = a[j];
        (*qs[group_index]) << partitions[j] << '\n';
    }
}

pllresult Optimiser::get_parameters(PLLUPtr&& pll) {
    int cap = pll->get_number_of_partitions();
    std::vector<perlocus> res(cap);
    for (int i=0; i < cap; ++i) {
        res[i].alpha = pll->get_alpha(i);
        res[i].freqs = pll->get_frequencies(i);
        res[i].rates = pll->get_rates(i);
        res[i].likelihood = (*pll)[i]->partitionLH;
    }
    return pllresult{res, pll->get_tree(), pll->get_likelihood()};
};
