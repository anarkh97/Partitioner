#include <mpi.h>
#include <ctime>
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <argparse/argparse.hpp>

#include "CommandLineData.hpp"
#include "LogUtils.h"

MPI_Comm comm;
double start_time;

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    start_time = LogUtils::WallTime();
    comm = MPI_COMM_WORLD;

    CommandLineData(comm, argc, argv);

    MPI_Finalize();
    return 0;

}
