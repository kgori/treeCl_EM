#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <pll/pll.h>
#include "memory_management.h"
#include "PLL.h"

using namespace std;
using namespace std::chrono;

string MYFILE="data/conc.phy";
string MYPART="data/conc.partitions.txt";

bool is_file(string filename) {
    ifstream fl(filename.c_str());
    bool result = true;
    if (!fl) {
        result = false;
    }
    fl.close();
    return result;
}

alignmentUPtr parse_alignment_file(string path) {
    if (!is_file(path)) {
        cerr << "Couldn't find the alignment file " << path << endl;
        throw exception();
    }
    alignmentUPtr alignment;
    AlignmentDeleter del;
    alignment = alignmentUPtr(pllParseAlignmentFile(PLL_FORMAT_PHYLIP, path.c_str()), del);
    if (!alignment) {
        //cout << "Trying to parse as fasta" << endl;
        alignment = alignmentUPtr(pllParseAlignmentFile(PLL_FORMAT_FASTA, path.c_str()), del);
    }
    if (!alignment) {
        cerr << "Couldn't parse the alignment at " << path << endl;
        throw exception();
    }
    return alignment;
}

queueUPtr parse_partitions(string partitions) {
    if (is_file(partitions)) {
        return queueUPtr(pllPartitionParse(partitions.c_str()), QueueDeleter());
    }
    else {
        return queueUPtr(pllPartitionParseString(partitions.c_str()), QueueDeleter());
    }
}

vector<string> readlines(string& filename) {
    ifstream infile(filename);
    vector<string> result;
    string line;
    if (!infile) {
        cerr << "Error opening file" << endl;
        return result;
    }

    while (std::getline(infile, line)) result.push_back(line);
    return result;
}

template<typename Iter>
string stringjoin(Iter it, Iter end, char delim = '\n') {
    ostringstream oss;
    for (;it != end - 1; ++it) {
        oss << *it << delim;
    }
    oss << *(it);
    return oss.str();
}

template<typename Iter>
void print_container(Iter it, Iter end, char delim = ' ') {
    for (;it != end - 1; ++it) {
        cout << *it << " ";
    }
    cout << *(it) << endl;
}

string get_tree(pllInstance* tr, partitionList* partitions) {
    pllTreeToNewick(tr->tree_string, tr, partitions,
                    tr->start->back, PLL_TRUE, PLL_TRUE,
                    PLL_FALSE, PLL_FALSE, PLL_FALSE,
                    0, PLL_FALSE, PLL_FALSE);
    string tree(tr->tree_string);
    tree.erase(std::remove(tree.begin(), tree.end(), '\n'), tree.end());
    return tree;
}

vector<int> get_within_group_index(const vector<int>& assignment) {
    auto max_elem = std::max_element(assignment.begin(), assignment.end());
    int ngrp = 1 + assignment[std::distance(assignment.begin(), max_elem)];
    vector<int> grpix(ngrp, 0);
    vector<int> result(assignment.size(), 0);
    for (int i = 0; i < assignment.size(); ++i) {
        int grp = assignment[i];
        result[i] = grpix[grp]++;
    }
    return result;
}


// Using mains for quick testing ATM
/*
int main() {
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    pllInstanceAttr attr;
    attr.rateHetModel = PLL_GAMMA;
    attr.fastScaling = PLL_FALSE;
    attr.saveMemory = PLL_FALSE;
    attr.useRecom = PLL_FALSE;
    attr.randomNumberSeed = 12345;
    attr.numberOfThreads = 8;

    vector<string> partitions = readlines(MYPART);
    queueUPtr q;
    alignmentUPtr al;

    for (int i = 0; i < partitions.size(); ++i) {
        al = parse_alignment_file(MYFILE);
        q = parse_partitions(partitions[i].c_str());
        auto pll = make_unique<PLL>(attr, q.get(), al.get());
        pll->tree_search(true);
        cout << "Likelihood = " << pll->get_likelihood() << endl;
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << duration;
    return 0;
}

*/
int main() {
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    pllInstanceAttr attr;
    attr.rateHetModel = PLL_GAMMA;
    attr.fastScaling = PLL_FALSE;
    attr.saveMemory = PLL_FALSE;
    attr.useRecom = PLL_FALSE;
    attr.randomNumberSeed = 12345;
    attr.numberOfThreads = 1;

    vector<string> partitions = readlines(MYPART);
    queueUPtr q;
    alignmentUPtr al;
    vector<PLLUPtr> insts;
    std::vector<std::thread> threads;

    for (int i = 0; i < partitions.size(); ++i) {
        al = parse_alignment_file(MYFILE);
        q = parse_partitions(partitions[i].c_str());
        insts.push_back(make_unique<PLL>(attr, q.get(), al.get()));
    }

    for (auto& pll : insts) {
        threads.push_back(pll->optimise_in_thread());
    }

    for (auto& t : threads) {
        t.join();
    }

    double lnlsum = 0;
    for (auto& pll : insts) {
        lnlsum += pll->get_likelihood();
        cout << pll->get_likelihood() << " " << pll->get_likelihood() / (*pll)[0]->width << endl;
    }
    cout << "LNLsum = " << lnlsum << endl;

    // Put in own scope - this could be a class / collection of functions
    {
        // BREAK INTO FN1 - Build vector of PLLs for each cluster
        // --------------
        // Assignments of loci to clusters (position in vector corresponds to locus,
        //                                  value at position corresponds to group index)
        vector<int> assgn{0, 0, 0, 1, 0, 2, 1, 2, 1, 1, 0, 2, 2, 0, 2};

        // Bookkeeping: within group index
        auto wgi = get_within_group_index(assgn);

        // Number of groups is 1 higher than maximum group index (assume 0-based indexing)
        auto max_elem = std::max_element(assgn.begin(), assgn.end());
        int numgrps = 1 + assgn[std::distance(assgn.begin(), max_elem)];

        // Initialise vector to hold PLL for each group
        vector<PLLUPtr> grps;
        grps.reserve(numgrps);

        // Initialise vector to build up partition info for each group
        vector<vector<string>> qs(numgrps);

        for (int i = 0; i < assgn.size(); ++i) {
            int group_index = assgn[i];
            qs[group_index].push_back(partitions[i]);
        }

        for (int j = 0; j < numgrps; ++j) {
            al = parse_alignment_file(MYFILE);
            string qstring = stringjoin(qs[j].begin(), qs[j].end()); // join together partition info from prev step
            q = parse_partitions(qstring.c_str());
            grps.push_back(make_unique<PLL>(attr, q.get(), al.get()));
        }
        // END FN1

        // BREAK INTO FN2 - Do computation with vector of PLLs
        // Initialise vector of threads to do PLL computation
        std::vector<std::thread> thrds;
        for (auto& pll : grps) {
            thrds.push_back(pll->optimise_in_thread()); // spawn thread
        }

        for (auto& t : thrds) {
            t.join(); // collect results
        }
        // END FN2

        // BREAK INTO FN3 - Process results
        // Process results (just printing for now)
        lnlsum = 0;
        for (auto& pll : grps) {
            lnlsum += pll->get_likelihood();
            cout << pll->get_likelihood() << endl;
        }
        cout << "LNLsum = " << lnlsum << endl;

        // These arrays (just printed ATM) will form basis for making moves between clusters
        for (auto& locus : insts) {
            for (auto& pll : grps) {
                locus->set_tree(pll->get_tree());
                cout << locus->get_likelihood() << " ";
            }
            cout << endl;
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
        vector<int> assgn{0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2};

        // Bookkeeping: within group index
        auto wgi = get_within_group_index(assgn);

        // Number of groups is 1 higher than maximum group index (assume 0-based indexing)
        auto max_elem = std::max_element(assgn.begin(), assgn.end());
        int numgrps = 1 + assgn[std::distance(assgn.begin(), max_elem)];

        // Initialise vector to hold PLL for each group
        vector<PLLUPtr> grps;
        grps.reserve(numgrps);

        // Initialise vector to build up partition info for each group
        vector<vector<string>> qs(numgrps);

        for (int i = 0; i < assgn.size(); ++i) {
            int group_index = assgn[i];
            qs[group_index].push_back(partitions[i]);
        }

        for (int j = 0; j < numgrps; ++j) {
            al = parse_alignment_file(MYFILE);
            string qstring = stringjoin(qs[j].begin(), qs[j].end()); // join together partition info from prev step
            q = parse_partitions(qstring.c_str());
            grps.push_back(make_unique<PLL>(attr, q.get(), al.get()));
        }
        // END FN1

        // BREAK INTO FN2 - Do computation with vector of PLLs
        // Initialise vector of threads to do PLL computation
        std::vector<std::thread> thrds;
        for (auto& pll : grps) {
            thrds.push_back(pll->optimise_in_thread()); // spawn thread
        }

        for (auto& t : thrds) {
            t.join(); // collect results
        }
        // END FN2

        // BREAK INTO FN3 - Process results
        // Process results (just printing for now)
        lnlsum = 0;
        for (auto& pll : grps) {
            lnlsum += pll->get_likelihood();
            cout << pll->get_likelihood() << endl;
        }
        cout << "LNLsum = " << lnlsum << endl;

        // These arrays (just printed ATM) will form basis for making moves between clusters
        for (auto& locus : insts) {
            for (auto& pll : grps) {
                locus->set_tree(pll->get_tree());
                cout << locus->get_likelihood() << " ";
            }
            cout << endl;
        }
        // END FN3

        // FN4 - Do reassignments
        // ...
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    cout << "Program runtime: " << duration / 1000.0 << "s";
    return 0;
}
