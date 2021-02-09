cmake_minimum_required(VERSION 3.14)
project(OPPlab1_upd)

set(CMAKE_CXX_STANDARD 14)

add_executable(OPPlab1_upd ver1.cpp test1.cpp test1.h      ver2.cpp test2.cpp test2.h)
#add_executable(OPPlab1_upd ver2.cpp ver2.h)
# Require MPI for this project:

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(OPPlab1_upd ${MPI_LIBRARIES})