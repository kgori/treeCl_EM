//
// Created by Kevin Gori on 02/06/15.
//

#include "threadpool.h"
thread_local work_stealing_queue* work_stealing_thread_pool::local_work_queue;
thread_local unsigned work_stealing_thread_pool::my_index;
thread_local std::unique_ptr<local_queue_type> local_queue_thread_pool;
