#ifndef __SPACE_DISCRETIZATION_H__
#define __SPACE_DISCRETIZATION_H__

#include <mpi.h>
#include <petscdmplex.h>
#include <petscpartitioner.h>
#include <petscsection.h>

#include "CommandLineData.h"
#include "GlobalTypeDefs.h"
#include "LogUtils.h"

class SpaceDiscretization
{
    MPI_Comm& comm;
    CommandLineData& cli_data;

    //!< DMPlex object
    DM dm;

    //!< Ranges (ranks)
    IntScalarType rank_point_start;  //!< DAG mesh points
    IntScalarType rank_point_end;

    IntScalarType rank_cell_start;
    IntScalarType rank_cell_end;

    IntScalarType rank_vertex_start;
    IntScalarType rank_vertex_end;

    IntScalarType rank_face_start;
    IntScalarType rank_face_end;

    IntScalarType rank_edge_start;
    IntScalarType rank_edge_end;

    //!< Ranges (patches)


    //!< Element patch layout
    IntVectorType elem2patch;

   public:
    //!< constructor
    explicit SpaceDiscretization(MPI_Comm& comm, CommandLineData& cli_data);
    //!< destructor (always use Destroy)
    ~SpaceDiscretization( ) {};

    //!< Methods
    void CreateMeshPatches( );

    //!< Clean-up
    void Destroy( );

    //!< Accesors
    inline void GetInternalCellIndexBounds(IntScalarType& cell_start,
                                           IntScalarType& cell_end) const
    {
        cell_start = rank_cell_start;
        cell_end = rank_cell_end;
    }

    inline void GetInternalFaceIndexBounds(IntScalarType& face_start,
                                           IntScalarType& face_end) const
    {
        face_start = rank_face_start;
        face_end = rank_face_end;
    }

    inline void GetInternalEdgeIndexBounds(IntScalarType& edge_start,
                                           IntScalarType& edge_end) const
    {
        edge_start = rank_edge_start;
        edge_end = rank_edge_end;
    }

    inline void GetInternalVertexIndexBounds(IntScalarType& vert_start,
                                             IntScalarType& vert_end) const
    {
        vert_start = rank_vertex_start;
        vert_end = rank_vertex_end;
    }

    //!< Output
    void OutputMesh( );
};

#endif
