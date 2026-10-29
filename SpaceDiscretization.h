#ifndef __SPACE_DISCRETIZATION_H__
#define __SPACE_DISCRETIZATION_H__

#include <mpi.h>
#include <petscdmplex.h>
#include <petscpartitioner.h>
#include <petscsection.h>
#include <vector>

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
    int_t rank_point_start;  //!< DAG mesh points
    int_t rank_point_end;

    int_t rank_cell_start;
    int_t rank_cell_end;

    int_t rank_vertex_start;
    int_t rank_vertex_end;

    int_t rank_face_start;
    int_t rank_face_end;

    int_t rank_edge_start;
    int_t rank_edge_end;

    //!< Ranges (patches)


    //!< Element patch layout
    std::vector<int_t> elem2patch;

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
    inline void GetInternalCellIndexBounds(int_t& cell_start,
                                           int_t& cell_end) const
    {
        cell_start = rank_cell_start;
        cell_end = rank_cell_end;
    }

    inline void GetInternalFaceIndexBounds(int_t& face_start,
                                           int_t& face_end) const
    {
        face_start = rank_face_start;
        face_end = rank_face_end;
    }

    inline void GetInternalEdgeIndexBounds(int_t& edge_start,
                                           int_t& edge_end) const
    {
        edge_start = rank_edge_start;
        edge_end = rank_edge_end;
    }

    inline void GetInternalVertexIndexBounds(int_t& vert_start,
                                             int_t& vert_end) const
    {
        vert_start = rank_vertex_start;
        vert_end = rank_vertex_end;
    }

    //!< Output
    void OutputMesh( );
};

#endif
