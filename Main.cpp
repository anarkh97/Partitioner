#include <mpi.h>
#include <ctime>
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <argparse/argparse.hpp>
#include <petscdmplex.h>
#include <petscsys.h>

#include "GlobalTypeDefs.h"
#include "LogUtils.h"
#include "CommandLineData.h"
#include "SpaceDiscretization.h"

MPI_Comm comm;
real_t start_time;

int main(int argc, char *argv[])
{

    //!<-----------------------------------------------------------------------
    //!< Initialize MPI
    //!<-----------------------------------------------------------------------
    MPI_Init(&argc, &argv);
    start_time = LogUtils::WallTime();
    comm = MPI_COMM_WORLD;

    //!<-----------------------------------------------------------------------
    //!< Parse the command line inputs (all ranks)
    //!<-----------------------------------------------------------------------
    CommandLineData cli_data(comm, argc, argv);

    //!<-----------------------------------------------------------------------
    //!< Initialize PETSc (only forward the program name)
    //!<-----------------------------------------------------------------------
    PETSC_COMM_WORLD = comm;
    int_t petsc_argc = 1;
    PetscInitialize(&petsc_argc, &argv, NULL, NULL);

    //!<-----------------------------------------------------------------------
    //!< Create discrete space
    //!<-----------------------------------------------------------------------
    SpaceDiscretization spd(comm, cli_data);
    spd.OutputMesh();


    //!<-----------------------------------------------------------------------
    //!< Finalize
    //!<-----------------------------------------------------------------------
    LogUtils::Print("\n");
    LogUtils::PrintSuccess("==========================");
    LogUtils::PrintSuccess("    NORMAL TERMINATION    "); 
    LogUtils::PrintSuccess("==========================");
    LogUtils::Print("\n");

    spd.Destroy();

    PetscFinalize();
    MPI_Finalize();

    return 0;

}
