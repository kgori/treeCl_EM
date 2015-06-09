#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
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

/*
 * Function wrapper for multithreaded execution in threadpool.
 *
 * Pass in a reference to a member function [that returns void], the number of threads, a vector of unique pointers to
 * the data items to operate on, and the arguments that will be forwarded to the function
 */
template <typename FunctionType, typename ValueType, typename ...Args>
void threadpool(FunctionType&& fn, unsigned nthreads, std::vector<std::unique_ptr<ValueType>>& data, Args&&... args) {
    work_stealing_thread_pool pool(nthreads);
    std::vector<std::future<void>> futures;

    for (auto& data_item : data) {
        auto bound_fn = std::bind(fn, data_item.get(), std::forward<Args>(args)...);
        futures.push_back(pool.submit(bound_fn));
    }

    for (auto& fut : futures) {
        fut.get();
    }
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
    return 1 + assignment[std::distance(assignment.begin(), max_elem)];
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

double reassign(std::vector<int> assignment,
                const std::vector<std::string>& partitions,
                std::vector<PLLUPtr>& insts,  // dropped const qualifier to work w/threadpool function
                std::vector<int>& new_assignment)
{
    int numgrps = get_number_of_groups(assignment);
    auto wgi = get_within_group_index(assignment, numgrps);
    std::vector<PLLUPtr> grps;
    grps.reserve(numgrps);
    std::vector<std::vector<std::string>> qs(numgrps);
    for (int i = 0; i < assignment.size(); ++i) {
        int group_index = assignment[i];
        qs[group_index].push_back(partitions[i]);
    }

    alignmentUPtr al;
    queueUPtr q;
    pllInstanceAttr attr = get_attr(PLL_GAMMA, PLL_FALSE, PLL_FALSE, PLL_FALSE, 12345);
    for (int j = 0; j < numgrps; ++j) {
        al = parse_alignment_file(MYFILE);
        std::string qstring = stringjoin(qs[j].begin(), qs[j].end()); // join together partition info from prev step
        q = parse_partitions(qstring.c_str());
        grps.push_back(std::make_unique<PLL>(attr, q.get(), al.get()));
    }

    //
    for (int i = 0; i < insts.size(); ++i) {
        grps[assignment[i]]->set_alpha(insts[i]->get_alpha(0), wgi[i], false);
        grps[assignment[i]]->set_frequencies(insts[i]->get_frequencies(0), wgi[i], false);
    }
    //

//    threadpool(&PLL::optimise, 1, grps, false, false, false, true, 0.1, false);
    threadpool(&PLL::tree_search, grps.size(), grps, false);
    double lnlsum = 0;
    std::vector<std::vector<double>> reassignment_matrix(insts.size(), std::vector<double>(numgrps, 0));

    for (int i=0; i<numgrps; ++i) {
        lnlsum += grps[i]->get_likelihood();
        std::string tree = grps[i]->get_tree();
        threadpool(&PLL::set_tree, 1, insts, tree); // SIGSEGV??
        for (int j=0; j < insts.size(); ++j) {
            reassignment_matrix[j][i] = insts[j]->get_likelihood();
        }
    }
    for (int i=0; i<insts.size(); ++i) {
        auto max_iter = std::max_element(reassignment_matrix[i].begin(), reassignment_matrix[i].end());
        new_assignment[i] = std::distance(reassignment_matrix[i].begin(), max_iter);
    }
    restricted_growth_notation(new_assignment);
    return lnlsum;
}

double strain(PLL* pll) {
    return pll->get_likelihood() / pll->sites();
}


// Store some results / observations here
struct parameters {
    double alpha;
    std::vector<double> freqs;
    std::vector<double> rates;
    double likelihood;
    double nsites;
    std::string tree;
};

int main()
{
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    pllInstanceAttr attr;
    attr.rateHetModel = PLL_GAMMA;
    attr.fastScaling = PLL_FALSE;
    attr.saveMemory = PLL_FALSE;
    attr.useRecom = PLL_FALSE;
    attr.randomNumberSeed = 12345;
    attr.numberOfThreads = std::thread::hardware_concurrency();


    std::vector<std::string> partitions = readlines(MYPART);
    unsigned nloci = partitions.size();
    queueUPtr q; // Can reuse these pointers
    alignmentUPtr al;
    std::vector<parameters> params;
    std::vector<PLLUPtr> insts;

    for (int i = 0; i < nloci; ++i) {
        al = parse_alignment_file(MYFILE);
        q = parse_partitions(partitions[i].c_str());
        PLLUPtr inst = std::make_unique<PLL>(attr, q.get(), al.get());
        inst->tree_search(true);
        parameters p;
        p.alpha = inst->get_alpha(0);
        p.freqs = inst->get_frequencies(0);
        p.rates = inst->get_rates(0);
        p.likelihood = inst->get_likelihood();
        p.nsites = (*inst)[0]->width;
        p.tree = inst->get_tree();
        params.push_back(p);
    }

    for (int i = 0; i < params.size(); ++i) {
        std::cout << i << ":" << params[i].alpha << std::endl;
    }

//    threadpool(&PLL::optimise, insts.size(), insts, true, true, true, true, 0.01, false);
//    threadpool(&PLL::tree_search, insts.size(), insts, true);

    std::vector<std::string> trees;

    //std::vector<double> dm;
    //for (int i = 0; i < insts.size(); ++i) {
        //for (int j = 0; j < insts.size(); ++j) {
            //insts[i]->set_tree(trees[j]);
            //dm.push_back(strain(insts[i].get()));
            //std::cout << strain(insts[i].get()) << " ";
        //}
        //std::cout << std::endl;
    //}

    //for (int i = 0; i < insts.size(); ++i) {
//        double refval = dm[i*insts.size()+i];
//        for (int j = 0; j < insts.size(); ++j) {
//            dm[i*insts.size()+j] -= refval;
//            std::cout << dm[i*insts.size()+j] << " ";
//        }
//        std::cout << std::endl;
//    }
//
//
//    std::vector<double> alphas;
//    double lnlsum = 0;
//    for (auto& pll : insts)
//    {
//        lnlsum += pll->get_likelihood();
//        std::cout << pll->get_likelihood() << " " << pll->get_likelihood() / (*pll)[0]->width << std::endl;
//        alphas.push_back(pll->get_alpha(0));
//    }
//    std::cout << "LNLsum = " << lnlsum << std::endl;
//    print_container(alphas.begin(), alphas.end());
//
//    std::vector<int> assgn{0, 1, 2, 3, 4, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
//    std::vector<int> new_assgn(assgn.size(), 0);
//    double bestlnl = PLL_UNLIKELY;
//    double newlnl = reassign(assgn, partitions, insts, new_assgn);
//    while(newlnl > bestlnl)
//    {
//        bestlnl = newlnl;
//        assgn = new_assgn;
//        newlnl = reassign(assgn, partitions, insts, new_assgn);
//        std::cout << "EOL:bestlnl = " << std::fixed << std::setw(11) << std::setprecision(5) << bestlnl << std::endl;
//        std::cout << "EOL:newlnl = " << std::fixed << std::setw(11) << std::setprecision(5) << newlnl << std::endl;
//    }
//    std::cout << stringjoin(assgn.begin(), assgn.end(), ' ') << std::endl;
//
//    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
//    std::cout << "Program runtime: " << duration / 1000.0 << "s" << std::endl;
//    return 0;
}
