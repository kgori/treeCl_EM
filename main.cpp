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
        threads.push_back(pll->tree_search_in_thread());
    }

    for (auto& t : threads) {
        t.join();
    }

    double lnlsum;
    for (auto& pll : insts) {
        lnlsum += pll->get_likelihood();
        cout << pll->get_likelihood() << " " << pll->get_likelihood() / (*pll)[0]->width << endl;
    }
    cout << "LNLsum = " << lnlsum << endl;

    // Put in own scope
    {
        // Assignments of loci to clusters (position in vector corresponds to locus,
        //                                  value at position corresponds to group index)
        vector<int> assgn{0, 0, 0, 0, 1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 0};

        // Number of groups is 1 higher than maximum group index (assume 0-based indexing)
        int numgrps = 1 + *(std::max_element(assgn.begin(), assgn.end()));
        cout << "Number of groups = " << numgrps << endl;

        // Initialise vector to hold PLL for each group
        vector<PLLUPtr> grps;
        grps.reserve(numgrps);
        cout << "Init & reserve grps" << endl;

        // Initialise vector to build up partition info for each group
        vector<vector<string>> qs(numgrps, vector<string>());
        cout << "Init qs to size " << qs.size() << endl;

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

        // Initialise vector of threads to do PLL computation
        std::vector<std::thread> thrds;
        for (auto& pll : grps) {
            thrds.push_back(pll->optimise_in_thread()); // spawn thread
        }

        for (auto& t : thrds) {
            t.join(); // collect results
        }

        // Process results (just printing for now)
        double lnlsum = 0;
        for (auto& pll : grps) {
            lnlsum += pll->get_likelihood();
            cout << pll->get_likelihood() << endl;
        }
        cout << "LNLsum = " << lnlsum << endl;
    }

    // Put in own scope
    {
        // Assignments of loci to clusters (position in vector corresponds to locus,
        //                                  value at position corresponds to group index)
        vector<int> assgn{0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2};

        // Number of groups is 1 higher than maximum group index (assume 0-based indexing)
        int numgrps = 1 + *(std::max_element(assgn.begin(), assgn.end()));
        cout << "Number of groups = " << numgrps << endl;

        // Initialise vector to hold PLL for each group
        vector<PLLUPtr> grps;
        grps.reserve(numgrps);
        cout << "Init & reserve grps" << endl;

        // Initialise vector to build up partition info for each group
        vector<vector<string>> qs(numgrps, vector<string>());
        cout << "Init qs to size " << qs.size() << endl;

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

        // Initialise vector of threads to do PLL computation
        std::vector<std::thread> thrds;
        for (auto& pll : grps) {
            thrds.push_back(pll->optimise_in_thread()); // spawn thread
        }

        for (auto& t : thrds) {
            t.join(); // collect results
        }

        // Process results (just printing for now)
        double lnlsum = 0;
        for (auto& pll : grps) {
            lnlsum += pll->get_likelihood();
            cout << pll->get_likelihood() << endl;
        }
        cout << "LNLsum = " << lnlsum << endl;
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
    cout << duration;
    return 0;
}


/*
int main() {
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
    q = parse_partitions(stringjoin(partitions.begin(), partitions.end(), '\n').c_str());
    al = parse_alignment_file(MYFILE);
    auto pll = make_unique<PLL>(attr, q.get(), al.get());
    std::vector<std::vector<int>> origWeights;
    for (int i = 0; i<pll->get_number_of_partitions(); ++i) {
        auto part = (*pll)[i];
        std::vector<int> v;
        v.reserve(part->width);
        v.assign(part->wgt, part->wgt + part->width);
        origWeights.push_back(v);
        print_container(v.begin(), v.end());
        cout << endl;
    }

    cout << "Likelihood = " << pll->get_likelihood() << endl;

    for (int i = 1; i < pll->get_number_of_partitions(); ++i) {
        auto part = (*pll)[i];
        for (int j = 0; j < part->width; ++j) {
            part->wgt[j] = 0;
        }
    }

    cout << "Likelihood = " << pll->get_likelihood() << endl;
    pll->tree_search(true);
    cout << "Likelihood = " << pll->get_likelihood() << endl;

    std::vector<std::string> trees;
    std::vector<PLLUPtr> insts;


    return 0;
    for (string& part : partitions) {
        al = parse_alignment_file(MYFILE);
//        q = parse_partitions(stringjoin(partitions).c_str());
        q = parse_partitions(part.c_str());
//        auto pll = make_unique<PLL>(attr, q.get(), al.get());
//        auto pll2 = std::move(pll);
        insts.push_back(make_unique<PLL>(attr, q.get(), al.get()));
    }

//    for (auto it = insts.begin(); it != insts.end(); ++it ) {
//        auto pll = it->get();
//        pll->tree_search(true);
//        cout << pll->get_likelihood() << " " << pll->get_likelihood() / (*pll)[0]->width << endl;
////        trees.push_back(pll->get_tree());
//    }

    std::vector<std::thread> threads;
    for (auto& pll : insts) {
        threads.push_back(pll->tree_search_in_thread());
    }

    for (auto& t : threads) {
        t.join();
    }

    for (auto& pll : insts) {
        cout << pll->get_likelihood() << " " << pll->get_likelihood() / (*pll)[0]->width << endl;
    }

//    q = parse_partitions(partitions[0].c_str());
//    auto pll = make_unique<PLL>(attr, q.get(), al.get());
//
//    pll->optimise(true, true, true, true, 0.1);
//    cout << pll->get_likelihood() << endl;
//    cout << pll->get_tree() << endl;
//    for (int i=0; i < pll->get_number_of_partitions(); ++i) {
//        cout << (*pll)[i]->width << " "
//             << (*pll)[i]->partitionLH << " "
//             << (*pll)[i]->partitionLH / (*pll)[i]->width << endl;
//    }

    return 0;
}
*/