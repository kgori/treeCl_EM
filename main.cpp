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

double strain(PLL* pll) {
    return pll->get_likelihood() / pll->sites();
}


std::mt19937 engine(12345);

int main()
{


    auto attr = std::make_shared<pllInstanceAttr>();
    attr->rateHetModel = PLL_GAMMA;
    attr->fastScaling = PLL_FALSE;
    attr->saveMemory = PLL_FALSE;
    attr->useRecom = PLL_FALSE;
    attr->randomNumberSeed = 12345;
    attr->numberOfThreads = std::thread::hardware_concurrency();

    std::vector<std::string> partitions = utils::readlines(MYPART);
    auto o = Optimiser(MYFILE, partitions, attr);
    o.set_assignment(std::vector<int>{0,1,2,0,1,2,0,1,2,0,1,2,0,1,2});
    std::vector<int> x = o.get_assignment();
    std::vector<double> y = o.get_proportions(1);
    utils::print_container(x.begin(), x.end());
    utils::print_container(y.begin(), y.end());

    o.mStep();
    o.eStep();
    std::cout << "VTAB" << std::endl;
    o.vtab->print();


    std::uniform_real_distribution<double> dist(0, 1);

    std::vector<double> a = utils::csum(utils::scale_by_sum(std::vector<double>{0.5, 1, 2}));
    utils::print_container(a.begin(), a.end());
    double u;
    size_t selection;
    std::vector<size_t> selections;
    for (int j=0; j<100000; ++j) {
        selection = utils::random_select(a);
        selections.push_back(selection);
    }


    for (int k=0; k < a.size()+1; ++k)
        std::cout << k << " chosen " << std::count(selections.begin(), selections.end(), k) << " times;" << std::endl;

    double vsmall = 9e-200;
    std::cout << pow(vsmall,0.01) << std::endl;
    std::cout << exp(0.01*log(vsmall)) << std::endl;
    return 0;
}
