cmake_minimum_required(VERSION 3.30)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

project(NPART LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-reorder")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

include_directories("${CMAKE_SOURCE_DIR}")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(MPI REQUIRED COMPONENTS CXX)

# fetch latest argparse
include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
)
FetchContent_MakeAvailable(json)

find_package(PETSc REQUIRED)
message(STATUS "Found PETSc: ${PETSc_LIBRARY} (version ${PETSc_VERSION})")
if(NOT "${PETSc_LIBRARY_DIR}" STREQUAL "")
  link_directories("${PETSc_LIBRARY_DIR}")
endif()
if(NOT "${PETSc_INCLUDE_DIR}" STREQUAL "")
  include_directories("${PETSc_INCLUDE_DIR}")
endif()

add_executable(npart
    Main.cpp
    LogUtils.cpp
    SpaceDiscretization.cpp
)

target_link_libraries(npart petsc mpi)
target_link_libraries(npart argparse)
target_link_libraries(npart nlohmann_json::nlohmann_json)
