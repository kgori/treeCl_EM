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
    }
};

struct QueueDeleter {
    void operator()(pllQueue *pQueue) {
        pllQueuePartitionsDestroy(&pQueue);
    }
};

struct NewickDeleter {
    void operator()(pllNewickTree *pNewick) {
        pllNewickParseDestroy(&pNewick);
    }
};

struct InstanceDeleter {
    void operator()(pllInstance *pInstance) {
        pllDestroyInstance(pInstance);
    }
};

using alignmentUPtr = std::unique_ptr<pllAlignmentData, AlignmentDeleter>;
using newickUPtr = std::unique_ptr<pllNewickTree, NewickDeleter>;
using queueUPtr = std::unique_ptr<pllQueue, QueueDeleter>;
using instanceUPtr = std::unique_ptr<pllInstance, InstanceDeleter>;
using attrSPtr = std::shared_ptr<pllInstanceAttr>;

#endif //TREECL_EM_MEMORY_MANAGEMENT_H
