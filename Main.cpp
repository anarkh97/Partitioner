#include <mpi.h>
#include <petscdmplex.h>
#include <petscsys.h>
#include <argparse/argparse.hpp>
#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include "CommandLineData.h"
#include "GlobalTypeDefs.h"
#include "LogUtils.h"
#include "SpaceDiscretization.h"

MPI_Comm comm;
RealScalarType start_time;

int main(int argc, char* argv[])
{
    //!<-----------------------------------------------------------------------
    //!< Initialize MPI
    //!<-----------------------------------------------------------------------
    MPI_Init(&argc, &argv);
    start_time = LogUtils::WallTime( );
    comm = MPI_COMM_WORLD;

    //!<-----------------------------------------------------------------------
    //!< Parse the command line inputs (all ranks)
    //!<-----------------------------------------------------------------------
    CommandLineData cli_data(comm, argc, argv);

    //!<-----------------------------------------------------------------------
    //!< Initialize PETSc (only forward the program name)
    //!<-----------------------------------------------------------------------
    PETSC_COMM_WORLD = comm;
    IntScalarType petsc_argc = 1;
    PetscInitialize(&petsc_argc, &argv, NULL, NULL);

    //!<-----------------------------------------------------------------------
    //!< Create discrete space
    //!<-----------------------------------------------------------------------
    SpaceDiscretization spd(comm, cli_data);
    spd.CreateMeshPatches( );
    spd.OutputMesh( );


    //!<-----------------------------------------------------------------------
    //!< Finalize
    //!<-----------------------------------------------------------------------
    LogUtils::Print("\n");
    LogUtils::PrintSuccess("==========================");
    LogUtils::PrintSuccess("    NORMAL TERMINATION    ");
    LogUtils::PrintSuccess("==========================");
    LogUtils::Print("\n");

    spd.Destroy( );

    PetscFinalize( );
    MPI_Finalize( );

    return 0;
}
