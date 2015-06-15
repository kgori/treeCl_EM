cmake_minimum_required(VERSION 3.2)
project(treeCl_EM)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -g")
set(MY_LIB_LINK_LIBRARIES -lpll-avx-pthreads)
add_subdirectory(data)

set(SOURCE_FILES
    memory_management.h
    PLL.cpp
    main.cpp)

add_executable(treeCl_EM ${SOURCE_FILES} Optimiser.cpp Optimiser.h utils.h utils.cpp ValueTable.cpp ValueTable.h)
TARGET_LINK_LIBRARIES(treeCl_EM ${MY_LIB_LINK_LIBRARIES})

