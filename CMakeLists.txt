cmake_minimum_required(VERSION 3.30)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

project(PARTITIONER C CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-reorder")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(MPI REQUIRED COMPONENTS CXX)

# fetch latest argparse
include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

include_directories("${CMAKE_SOURCE_DIR}")

add_executable(partitioner
    Main.cpp
    LogUtils.cpp
)

target_link_libraries(partitioner mpi)
target_link_libraries(partitioner argparse)
