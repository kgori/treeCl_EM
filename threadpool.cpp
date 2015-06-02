//
// Created by Kevin Gori on 02/06/15.
//

#include "threadpool.h"
thread_local work_stealing_queue* thread_pool::local_work_queue;
thread_local unsigned thread_pool::my_index;

//TEST FN

