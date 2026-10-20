#ifndef __SPACE_DISCRETIZATION_H__
#define __SPACE_DISCRETIZATION_H__

#include <mpi.h>
#include <vector>
#include <petscdmplex.h>
#include <petscsection.h>
#include <petscpartitioner.h>

#include "GlobalTypeDefs.h"
#include "LogUtils.h"
#include "CommandLineData.h"

class SpaceDiscretization {

    MPI_Comm& comm;
    CommandLineData& cli_data;

    //!< DMPlex object
    DM dm;

    //!< Hypergraph
    std::vector<int_t> elem_node_offsets;
    std::vector<int_t> elem_nodes;
    std::vector<byte_t> node_bytes;

    //!< Element adjacency
    std::vector<int_t> elem_adj_offsets;
    std::vector<int_t> elem_adj;
    std::vector<int_t> shared_node_counts;

public:

    //!< constructor
    explicit SpaceDiscretization(MPI_Comm& comm, CommandLineData& cli_data);
    //!< destructor (always use Destroy)
    ~SpaceDiscretization() { };

    //!< Clean-up
    void Destroy();
    
    //!< Output
    void OutputMesh();

};

#endif
