#include "SpaceDiscretization.h"
#include "CommandLineData.h"
#include "petscdm.h"

//-----------------------------------------------------------------------------

SpaceDiscretization::SpaceDiscretization(MPI_Comm& _comm,
                                         CommandLineData& _cli_data)
    : comm(_comm), cli_data(_cli_data)
{

    //!<-----------------------------------------------------------------------
    //!< Load mesh on rank 0 and then partition
    //!<-----------------------------------------------------------------------
    DM dm_dist; //!< Tempory DM to hold per-rank mesh
    
    DMPlexCreateGmshFromFile(comm, cli_data.mesh_file.c_str(), true, &dm);
    DMPlexDistribute(dm, 1, NULL, &dm_dist);
    DMDestroy(&dm);
    dm = dm_dist;

}

//-----------------------------------------------------------------------------

void
SpaceDiscretization::Destroy()
{
    DMDestroy(&dm);
}

//-----------------------------------------------------------------------------

void
SpaceDiscretization::OutputMesh()
{
    //!< Not implemented yet.
}

//-----------------------------------------------------------------------------
