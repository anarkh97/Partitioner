#include "SpaceDiscretization.h"
#include "CommandLineData.h"
#include "petscdm.h"
#include "petscdmplex.h"
#include "petscpartitioner.h"
#include "petscsection.h"
#include "petscsystypes.h"
#include "petscvec.h"

//-----------------------------------------------------------------------------

SpaceDiscretization::SpaceDiscretization(MPI_Comm& _comm,
                                         CommandLineData& _cli_data)
    : comm(_comm), cli_data(_cli_data)
{

    int_t size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    //!<-----------------------------------------------------------------------
    //!< Load mesh on rank 0 and then partition
    //!<-----------------------------------------------------------------------
    DMPlexCreateGmshFromFile(comm, cli_data.mesh_file.c_str(), true, &dm);

    if (size > 1) {
        DM dm_dist; //!< Tempory DM to hold per-rank mesh

        //!< Add cell weights (number of nodes per elements)
        int_t point_start;
        int_t point_end;

        int_t vertex_start;
        int_t vertex_end;
        int_t vertex_index;

        DMPlexGetChart(dm, &point_start, &point_end);
        DMPlexGetDepthStratum(dm, 0, &vertex_start, &vertex_end);

        PetscSection weights;
        PetscSectionCreate(comm, &weights);
        PetscSectionSetChart(weights, point_start, point_end);
        for (vertex_index = vertex_start; vertex_index < vertex_end; vertex_index++) {
            PetscSectionSetDof(weights, vertex_index, 1);
        }
        PetscSectionSetUp(weights);

        DMSetLocalSection(dm, weights);
        PetscSectionDestroy(&weights);

        PetscPartitioner partitioner;
        DMPlexGetPartitioner(dm, &partitioner);
        PetscPartitionerSetType(partitioner, PETSCPARTITIONERPARMETIS);
        PetscPartitionerSetFromOptions(partitioner);
        DMPlexDistribute(dm, 0, NULL, &dm_dist);
        DMDestroy(&dm);

        dm = dm_dist;
    }

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

    int_t rank;
    MPI_Comm_rank(comm, &rank);

    //!<-----------------------------------------------------------------------
    //!< Step 1: Create cell-wise vector
    //!<-----------------------------------------------------------------------
    int_t point_start;
    int_t point_end;

    int_t cell_start;
    int_t cell_end;
    int_t cell_index;

    DMPlexGetChart(dm, &point_start, &point_end);
    DMPlexGetHeightStratum(dm, 0, &cell_start, &cell_end);

    PetscSection section;
    PetscSectionCreate(comm, &section);
    PetscSectionSetNumFields(section, 1);
    PetscSectionSetFieldName(section, 0, "");
    PetscSectionSetChart(section, point_start, point_end);
    for (cell_index = cell_start; cell_index < cell_end; cell_index++) {
        PetscSectionSetDof(section, cell_index, 1);
        PetscSectionSetFieldDof(section, cell_index, 0, 1);
    }
    PetscSectionSetUp(section);

    DM dm_colors;
    DMClone(dm, &dm_colors);
    DMSetLocalSection(dm_colors, section);

    //!<-----------------------------------------------------------------------
    //!< Step 2: Create local vector
    //!<-----------------------------------------------------------------------
    real_t* colors_data;
    Vec     colors;

    DMCreateLocalVector(dm_colors, &colors);
    PetscObjectSetName((PetscObject)colors, "colors");

    //!<-----------------------------------------------------------------------
    //!< Step 3: Fill local vector with ranks
    //!<-----------------------------------------------------------------------
    VecGetArray(colors, &colors_data);

    for (cell_index = cell_start; cell_index < cell_end; cell_index++)
    {
        colors_data[cell_index-cell_start] = rank;
    }

    VecRestoreArray(colors, &colors_data);

    //!<-----------------------------------------------------------------------
    //!< Step 4: Create viewer and visualize mesh
    //!<-----------------------------------------------------------------------
    PetscViewer viewer;
    PetscViewerVTKOpen(comm, cli_data.out_file.c_str(), FILE_MODE_WRITE, &viewer);
    VecView(colors, viewer);

    //!<-----------------------------------------------------------------------
    //!< Step 4: Cleanup
    //!<-----------------------------------------------------------------------
    PetscSectionDestroy(&section);
    PetscViewerDestroy(&viewer);
    VecDestroy(&colors);
    DMDestroy(&dm_colors);

}

//-----------------------------------------------------------------------------
