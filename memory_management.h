//
// Created by Kevin Gori on 27/05/15.
//
#ifndef TREECL_EM_MEMORY_MANAGEMENT_H
#define TREECL_EM_MEMORY_MANAGEMENT_H

extern "C" {
#include "pll/pll.h"
}
#include <iostream>
#include <memory>

struct AlignmentDeleter {
    void operator()(pllAlignmentData *pAlignment) {
        pllAlignmentDataDestroy(pAlignment);
        //std::cout << "Alignment destroyed" << std::endl;
    }
};

struct QueueDeleter {
    void operator()(pllQueue *pQueue) {
        pllQueuePartitionsDestroy(&pQueue);
        //std::cout << "Queue destroyed" << std::endl;
    }
};

struct NewickDeleter {
    void operator()(pllNewickTree *pNewick) {
        pllNewickParseDestroy(&pNewick);
        //std::cout << "Newick destroyed" << std::endl;
    }
};

struct InstanceDeleter {
    void operator()(pllInstance *pInstance) {
        pllDestroyInstance(pInstance);
        //std::cout << "Instance Destroyed" << std::endl;
    }
};

typedef std::unique_ptr<pllAlignmentData, AlignmentDeleter> alignmentUPtr;
typedef std::unique_ptr<pllNewickTree, NewickDeleter> newickUPtr;
typedef std::unique_ptr<pllQueue, QueueDeleter> queueUPtr;
typedef std::unique_ptr<pllInstance, InstanceDeleter> instanceUPtr;

#endif //TREECL_EM_MEMORY_MANAGEMENT_H
