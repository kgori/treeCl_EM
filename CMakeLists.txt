cmake_minimum_required(VERSION 3.2)
project(treeCl_EM)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -g")
set(MY_LIB_LINK_LIBRARIES -lpll-avx)
add_subdirectory(data)

set(SOURCE_FILES
    PLL.cpp
    main.cpp)

add_executable(treeCl_EM ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(treeCl_EM ${MY_LIB_LINK_LIBRARIES})

