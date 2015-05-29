#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <pll/pll.h>
#include "memory_management.h"
#include "PLL.h"

using namespace std;

string MYFILE="/Users/kgori/Documents/repositories/treeCl_EM/data/conc.phy";
string MYPART="/Users/kgori/Documents/repositories/treeCl_EM/data/conc.partitions.txt";

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

string join(vector<string>& pieces, char delim='\n') {
    ostringstream oss;
    for (string& s : pieces) {
        oss << s << delim;
    }
    return oss.str();
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


//vector<string>

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
    q = parse_partitions(join(partitions).c_str());

    std::vector<std::string> trees;
    std::vector<PLLUPtr> insts;

    for (string& part : partitions) {
        auto al = parse_alignment_file(MYFILE);
//        q = parse_partitions(join(partitions).c_str());
        q = parse_partitions(part.c_str());
//        auto pll = make_unique<PLL>(attr, q.get(), al.get());
//        auto pll2 = std::move(pll);
        insts.push_back(make_unique<PLL>(attr, q.get(), al.get()));
    }

    for (auto it = insts.begin(); it != insts.end(); ++it ) {
        auto pll = it->get();
        pll->optimise(false, false, false, true);
        cout << pll->get_likelihood() << " " << pll->get_likelihood() / (*pll)[0]->width << endl;
        trees.push_back(pll->get_tree());
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
