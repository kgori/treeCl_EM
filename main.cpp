#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <pll/pll.h>
#include "PLL.h"


std::string MYFILE="data/conc.phy";
std::string MYPART="data/conc.partitions.txt";

bool is_file(std::string filename) {
    std::ifstream fl(filename.c_str());
    bool result = true;
    if (!fl) {
        result = false;
    }
    fl.close();
    return result;
}

alignmentUPtr parse_alignment_file(std::string path) {
    if (!is_file(path)) {
        std::cerr << "Couldn't find the alignment file " << path << std::endl;
        throw std::exception();
    }
    alignmentUPtr alignment;
    AlignmentDeleter del;
    alignment = alignmentUPtr(pllParseAlignmentFile(PLL_FORMAT_PHYLIP, path.c_str()), del);
    if (!alignment) {
        //std::cout << "Trying to parse as fasta" << std::endl;
        alignment = alignmentUPtr(pllParseAlignmentFile(PLL_FORMAT_FASTA, path.c_str()), del);
    }
    if (!alignment) {
        std::cerr << "Couldn't parse the alignment at " << path << std::endl;
        throw std::exception();
    }
    return alignment;
}

queueUPtr parse_partitions(std::string partitions) {
    if (is_file(partitions)) {
        return queueUPtr(pllPartitionParse(partitions.c_str()), QueueDeleter());
    }
    else {
        return queueUPtr(pllPartitionParseString(partitions.c_str()), QueueDeleter());
    }
}

std::vector<std::string> readlines(const std::string& filename) {
    std::ifstream infile(filename);
    std::vector<std::string> result;
    std::string line;
    if (!infile) {
        std::cerr << "Error opening file" << std::endl;
        return result;
    }

    while (std::getline(infile, line)) result.push_back(line);
    return result;
}

template<typename Iter>
std::string stringjoin(Iter it, Iter end, char delim = '\n') {
    std::ostringstream oss;
    if (it == end) return oss.str();
    --end;
    for (;it != end; ++it) {
        oss << *it << delim;
    }
    oss << *(it);
    return oss.str();
}

template<typename Iter>
void print_container(Iter it, Iter end, char delim = ' ') {
    if (it == end) return;
    --end;
    for (;it != end; ++it) {
        std::cout << *it << " ";
    }
    std::cout << *(it) << std::endl;
}

std::string get_tree(pllInstance* tr, partitionList* partitions) {
    pllTreeToNewick(tr->tree_string, tr, partitions,
                    tr->start->back, PLL_TRUE, PLL_TRUE,
                    PLL_FALSE, PLL_FALSE, PLL_FALSE,
                    0, PLL_FALSE, PLL_FALSE);
    std::string tree(tr->tree_string);
    tree.erase(std::remove(tree.begin(), tree.end(), '\n'), tree.end());
    return tree;
}

std::vector<int> get_within_group_index(const std::vector<int>& assignment) {
    auto max_elem = std::max_element(assignment.begin(), assignment.end());
    int ngrp = 1 + assignment[std::distance(assignment.begin(), max_elem)];
    std::vector<int> grpix(ngrp, 0);
    std::vector<int> result(assignment.size(), 0);
    for (int i = 0; i < assignment.size(); ++i) {
        int grp = assignment[i];
        result[i] = grpix[grp]++;
    }
    return result;
}

/*
 * Function wrapper for multithreaded execution in threadpool.
 *
 * Pass in a reference to a member function [that returns void], the number of threads, a vector of unique pointers to
 * the data items to operate on, and the arguments that will be forwarded to the function
 */
template <typename FunctionType, typename ValueType, typename ...Args>
void threadpool(FunctionType&& fn, unsigned nthreads, std::vector<std::unique_ptr<ValueType>>& data, Args&&... args) {
    thread_pool pool(nthreads);
    std::vector<std::future<void>> futures;

    for (auto& data_item : data) {
        auto bound_fn = std::bind(fn, data_item.get(), std::forward<Args>(args)...);
        futures.push_back(pool.submit(bound_fn));
    }

    for (auto& fut : futures) {
        fut.get();
    }
}

int main() {
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    pllInstanceAttr attr;
    attr.rateHetModel = PLL_GAMMA;
    attr.fastScaling = PLL_FALSE;
    attr.saveMemory = PLL_FALSE;
    attr.useRecom = PLL_FALSE;
    attr.randomNumberSeed = 12345;
    attr.numberOfThreads = 1;

    std::vector<std::string> partitions = readlines(MYPART);
    queueUPtr q;
    alignmentUPtr al;
    std::vector<PLLUPtr> insts;

    for (int i = 0; i < partitions.size(); ++i) {
        al = parse_alignment_file(MYFILE);
        q = parse_partitions(partitions[i].c_str());
        insts.push_back(std::make_unique<PLL>(attr, q.get(), al.get()));
    }

    threadpool(&PLL::optimise, 8, insts, true, true, true, true, 0.01, false);

    double lnlsum = 0;
    for (auto& pll : insts) {
        lnlsum += pll->get_likelihood();
        std::cout << pll->get_likelihood() << " " << pll->get_likelihood() / (*pll)[0]->width << std::endl;
    }
    std::cout << "LNLsum = " << lnlsum << std::endl;

    // Put in own scope - this could be a class / collection of functions
    {
        // BREAK INTO FN1 - Build vector of PLLs for each cluster
        // --------------
        // Assignments of loci to clusters (position in vector corresponds to locus,
        //                                  value at position corresponds to group index)
        std::vector<int> assgn{0, 0, 0, 1, 0, 2, 1, 2, 1, 1, 0, 2, 2, 0, 2};

        // Bookkeeping: within group index
        auto wgi = get_within_group_index(assgn);

        // Number of groups is 1 higher than maximum group index (assume 0-based indexing)
        auto max_elem = std::max_element(assgn.begin(), assgn.end());
        int numgrps = 1 + assgn[std::distance(assgn.begin(), max_elem)];

        // Initialise vector to hold PLL for each group
        std::vector<PLLUPtr> grps;
        grps.reserve(numgrps);

        // Initialise vector to build up partition info for each group
        std::vector<std::vector<std::string>> qs(numgrps);

        for (int i = 0; i < assgn.size(); ++i) {
            int group_index = assgn[i];
            qs[group_index].push_back(partitions[i]);
        }

        for (int j = 0; j < numgrps; ++j) {
            al = parse_alignment_file(MYFILE);
            std::string qstring = stringjoin(qs[j].begin(), qs[j].end()); // join together partition info from prev step
            q = parse_partitions(qstring.c_str());
            grps.push_back(std::make_unique<PLL>(attr, q.get(), al.get()));
        }
        // END FN1

        // BREAK INTO FN2 - Do computation with vector of PLLs
        // Initialise vector of threads to do PLL computation
        threadpool(&PLL::optimise, grps.size(), grps, false, false, false, true, 0.001, false);
        // END FN2

        // BREAK INTO FN3 - Process results
        // Process results (just printing for now)
        lnlsum = 0;
        for (auto& pll : grps) {
            lnlsum += pll->get_likelihood();
            std::cout << pll->get_likelihood() << std::endl;
        }
        std::cout << "LNLsum = " << lnlsum << std::endl;

        // These arrays (just printed ATM) will form basis for making moves between clusters
        for (auto& locus : insts) {
            for (auto& pll : grps) {
                locus->set_tree(pll->get_tree());
                std::cout << locus->get_likelihood() << " ";
            }
            std::cout << std::endl;
        }
        // END FN3

        // FN4 - Do reassignments
        // ...
    }

    // Put in own scope - this could be a class / collection of functions
    {
        // BREAK INTO FN1 - Build vector of PLLs for each cluster
        // --------------
        // Assignments of loci to clusters (position in vector corresponds to locus,
        //                                  value at position corresponds to group index)
        std::vector<int> assgn{0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2};

        // Bookkeeping: within group index
        auto wgi = get_within_group_index(assgn);

        // Number of groups is 1 higher than maximum group index (assume 0-based indexing)
        auto max_elem = std::max_element(assgn.begin(), assgn.end());
        int numgrps = 1 + assgn[std::distance(assgn.begin(), max_elem)];

        // Initialise vector to hold PLL for each group
        std::vector<PLLUPtr> grps;
        grps.reserve(numgrps);

        // Initialise vector to build up partition info for each group
        std::vector<std::vector<std::string>> qs(numgrps);

        for (int i = 0; i < assgn.size(); ++i) {
            int group_index = assgn[i];
            qs[group_index].push_back(partitions[i]);
        }

        for (int j = 0; j < numgrps; ++j) {
            al = parse_alignment_file(MYFILE);
            std::string qstring = stringjoin(qs[j].begin(), qs[j].end()); // join together partition info from prev step
            q = parse_partitions(qstring.c_str());
            grps.push_back(std::make_unique<PLL>(attr, q.get(), al.get()));
        }
        // END FN1

        // BREAK INTO FN2 - Do computation with vector of PLLs
        // Initialise vector of threads to do PLL computation
        threadpool(&PLL::optimise, grps.size(), grps, false, false, false, true, 0.001, false);
        // END FN2

        // BREAK INTO FN3 - Process results
        // Process results (just printing for now)
        lnlsum = 0;
        for (auto& pll : grps) {
            lnlsum += pll->get_likelihood();
            std::cout << pll->get_likelihood() << std::endl;
        }
        std::cout << "LNLsum = " << lnlsum << std::endl;

        // These arrays (just printed ATM) will form basis for making moves between clusters
        for (auto& locus : insts) {
            for (auto& pll : grps) {
                locus->set_tree(pll->get_tree());
                std::cout << locus->get_likelihood() << " ";
            }
            std::cout << std::endl;
        }
        // END FN3

        // FN4 - Do reassignments
        // ...
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    std::cout << "Program runtime: " << duration / 1000.0 << "s" << std::endl;
    return 0;
}
