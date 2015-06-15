#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <pll/pll.h>
#include <thread>
#include "PLL.h"
#include "Optimiser.h"
#include "utils.h"


std::string MYFILE="data/conc.phy";
std::string MYPART="data/conc.partitions.txt";




std::string get_tree(pllInstance* tr, partitionList* partitions) {
    pllTreeToNewick(tr->tree_string, tr, partitions,
                    tr->start->back, PLL_TRUE, PLL_TRUE,
                    PLL_FALSE, PLL_FALSE, PLL_FALSE,
                    0, PLL_FALSE, PLL_FALSE);
    std::string tree(tr->tree_string);
    tree.erase(std::remove(tree.begin(), tree.end(), '\n'), tree.end());
    return tree;
}



std::vector<int> get_within_group_index(const std::vector<int>& assignment, int numgroups) {
    std::vector<int> grpix(numgroups, 0);
    std::vector<int> result(assignment.size(), 0);
    for (int i = 0; i < assignment.size(); ++i) {
        int grp = assignment[i];
        result[i] = grpix[grp]++;
    }
    return result;
}
/*
 * Rewrite a partition vector in restricted growth notation (i.e. the first group is called '0',
 * the second is '1', ..., up to 'N-1'.
 */
void restricted_growth_notation(std::vector<int>& assignment) {
    int grpnum = 0;
    std::map<int, int> grpconv;
    for (auto it=assignment.begin(); it!=assignment.end(); it++) {
        auto search = grpconv.find(*it);
        if(search == grpconv.end()) {
            grpconv[*it] = grpnum;
            *it = grpnum++;
        }
        else {
            *it = search->second;
        }
    }
}

int get_number_of_groups(const std::vector<int>& assignment) {
    auto max_elem = std::max_element(assignment.begin(), assignment.end());
    return 1 + *max_elem;
}

pllInstanceAttr get_attr(int rateHetModel, int fastScaling, int saveMemory, int useRecom, int randomNumberSeed) {
    pllInstanceAttr attr;
    attr.rateHetModel = rateHetModel;
    attr.fastScaling = fastScaling;
    attr.saveMemory = saveMemory;
    attr.useRecom = useRecom;
    attr.randomNumberSeed = randomNumberSeed;
    attr.numberOfThreads = 1;
    return attr;
}

//double reassign(std::vector<int> assignment,
//                const std::vector<std::string>& partitions,
//                std::vector<PLLUPtr>& insts,  // dropped const qualifier to work w/threadpool function
//                std::vector<int>& new_assignment)
//{
//    int numgrps = get_number_of_groups(assignment);
//    auto wgi = get_within_group_index(assignment, numgrps);
//    std::vector<PLLUPtr> grps;
//    grps.reserve(numgrps);
//    std::vector<std::vector<std::string>> qs(numgrps);
//    for (int i = 0; i < assignment.size(); ++i) {
//        int group_index = assignment[i];
//        qs[group_index].push_back(partitions[i]);
//    }
//
//    alignmentUPtr al;
//    queueUPtr q;
//    pllInstanceAttr attr = get_attr(PLL_GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);
//    for (int j = 0; j < numgrps; ++j) {
//        al = parse_alignment_file(MYFILE);
//        std::string qstring = stringjoin(qs[j].begin(), qs[j].end()); // join together partition info from prev step
//        q = parse_partitions(qstring.c_str());
//        grps.push_back(std::make_unique<PLL>(attr, q.get(), al.get()));
//    }
//
//    //
//    for (int i = 0; i < insts.size(); ++i) {
//        grps[assignment[i]]->set_alpha(insts[i]->get_alpha(0), wgi[i], false);
//        grps[assignment[i]]->set_frequencies(insts[i]->get_frequencies(0), wgi[i], false);
//    }
//    //
//
////    threadpool(&PLL::optimise, 1, grps, false, false, false, true, 0.1, false);
//    //threadpool(&PLL::tree_search, grps.size(), grps, false);
//    double lnlsum = 0;
//    std::vector<std::vector<double>> reassignment_matrix(insts.size(), std::vector<double>(numgrps, 0));
//
//    for (int i=0; i<numgrps; ++i) {
//        lnlsum += grps[i]->get_likelihood();
//        std::string tree = grps[i]->get_tree();
//        //threadpool(&PLL::set_tree, 1, insts, tree); // SIGSEGV??
//        for (int j=0; j < insts.size(); ++j) {
//            reassignment_matrix[j][i] = insts[j]->get_likelihood();
//        }
//    }
//    for (int i=0; i<insts.size(); ++i) {
//        auto max_iter = std::max_element(reassignment_matrix[i].begin(), reassignment_matrix[i].end());
//        new_assignment[i] = std::distance(reassignment_matrix[i].begin(), max_iter);
//    }
//    restricted_growth_notation(new_assignment);
//    return lnlsum;
//}

double strain(PLL* pll) {
    return pll->get_likelihood() / pll->sites();
}


int main()
{
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    auto attr = std::make_shared<pllInstanceAttr>();
    attr->rateHetModel = PLL_GAMMA;
    attr->fastScaling = PLL_FALSE;
    attr->saveMemory = PLL_FALSE;
    attr->useRecom = PLL_FALSE;
    attr->randomNumberSeed = 12345;
    attr->numberOfThreads = std::thread::hardware_concurrency();

    std::vector<std::string> partitions = utils::readlines(MYPART);
    auto o = Optimiser(MYFILE, partitions, attr);
    o.set_assignment(std::vector<int>{0,0,0,0,0,2,2,2,2,2,1,1,1,1,1});
    std::vector<int> x = o.get_assignment();
    std::vector<double> y = o.get_proportions(1);
    utils::print_container(x.begin(), x.end());
    utils::print_container(y.begin(), y.end());
    o.mStep();
    std::cout << o.get_likelihood() << std::endl;
    return 0;
}
