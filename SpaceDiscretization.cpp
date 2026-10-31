#include "SpaceDiscretization.h"
#include "CommandLineData.h"
#include "GlobalTypeDefs.h"
#include "LogUtils.h"
#include "petscdm.h"
#include "petscdmplex.h"
#include "petscis.h"
#include "petscpartitioner.h"
#include "petscsection.h"
#include "petscsectiontypes.h"
#include "petscsys.h"
#include "petscsystypes.h"
#include "petscvec.h"

#include <algorithm>
#include <cmath>
#include <numeric>

using std::accumulate;
using std::ceil;
using std::iota;
using std::max;
using std::min;
using std::stable_sort;

constexpr IntScalarType INITIAL_SIZE = 1024;

//-----------------------------------------------------------------------------

SpaceDiscretization::SpaceDiscretization(MPI_Comm& _comm,
                                         CommandLineData& _cli_data)
    : comm(_comm), cli_data(_cli_data)
{
    Int32ScalarType size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    //!<-----------------------------------------------------------------------
    //!< Load mesh on rank 0 and then partition (rank-wise)
    //!<-----------------------------------------------------------------------
    DMPlexCreateGmshFromFile(comm, cli_data.mesh_file.c_str( ), true, &dm);

    IntScalarType dim;
    DMGetDimension(dm, &dim);

    if (dim != 3)
    {
        LogUtils::PrintError("Only 3D meshes are supported.");
        PetscFinalize( );
        MPI_Abort(comm, 0);
    }

    if (size > 1)
    {
        //!< Add cell weights (number of vertex per elements)
        IntScalarType point_start;
        IntScalarType point_end;

        IntScalarType vertex_start;
        IntScalarType vertex_end;
        IntScalarType vertex_index;

        DMPlexGetChart(dm, &point_start, &point_end);
        DMPlexGetDepthStratum(dm, 0, &vertex_start, &vertex_end);

        PetscSection weights;
        PetscSectionCreate(comm, &weights);
        PetscSectionSetChart(weights, point_start, point_end);
        for (vertex_index = vertex_start; vertex_index < vertex_end;
             vertex_index++)
        {
            PetscSectionSetDof(weights, vertex_index, 1);
        }
        PetscSectionSetUp(weights);

        DMSetLocalSection(dm, weights);
        PetscSectionDestroy(&weights);

        PetscPartitioner partitioner;
        DMPlexGetPartitioner(dm, &partitioner);
        PetscPartitionerSetType(partitioner, PETSCPARTITIONERPARMETIS);
        PetscPartitionerSetFromOptions(partitioner);

        DM dm_dist;  //!< Tempory DM to hold per-rank mesh
        DMPlexDistribute(dm, 0, NULL, &dm_dist);
        DMDestroy(&dm);

        //!< Update the DM object; now holds per-rank quantities
        dm = dm_dist;

        //!< Remove the weight section
        DMSetLocalSection(dm, NULL);
    }

    //!<-----------------------------------------------------------------------
    //!< Rank local ranges
    //!<-----------------------------------------------------------------------
    DMPlexGetChart(dm, &rank_point_start, &rank_point_end);
    DMPlexGetDepthStratum(dm, 0, &rank_vertex_start, &rank_vertex_end);
    DMPlexGetDepthStratum(dm, 1, &rank_edge_start, &rank_edge_end);
    DMPlexGetDepthStratum(dm, 2, &rank_face_start, &rank_face_end);
    DMPlexGetHeightStratum(dm, 0, &rank_cell_start, &rank_cell_end);
}

//-----------------------------------------------------------------------------

//!< AN: This can be something like, `AttachEquations` that takes in a vector
//!< of equations, gets the DOFs per stratum that an/all equations need, create
//!< patches based on required blocks, and then create `PetscSection` for each
//!< equation. We can then also create a global solution vector(s).
void SpaceDiscretization::CreateMeshPatches( )
{
    //!<-----------------------------------------------------------------------
    //!< Bytes required per DOF (hard-coded for now)
    //!< Currently for P2/Q2 elements; linear grid (xyz on vertex only).
    //!<-----------------------------------------------------------------------
    constexpr Int64ScalarType BYTE_PER_VERT =
        3 * sizeof(RealScalarType)     //!< Vector unknowns
        + 3 * sizeof(RealScalarType);  //!< Coordinates
    constexpr Int64ScalarType BYTE_PER_EDGE =
        3 * sizeof(RealScalarType);  //!< Vector unknowns
    constexpr Int64ScalarType BYTE_PER_FACE = 0;
    constexpr Int64ScalarType BYTE_PER_CELL = 0;

    //!< Helpers
    IntArrayType<4> stratum_starts = {rank_vertex_start, rank_edge_start,
                                      rank_face_start, rank_cell_start};
    IntArrayType<4> stratum_ends = {rank_vertex_end, rank_edge_end,
                                    rank_face_end, rank_cell_end};
    IntArrayType<4> stratum_size = {
        rank_vertex_end - rank_vertex_start, rank_edge_end - rank_edge_start,
        rank_face_end - rank_face_start, rank_cell_end - rank_cell_start};
    Int64ArrayType<4> bytes_per_stratum = {BYTE_PER_VERT, BYTE_PER_EDGE,
                                           BYTE_PER_FACE, BYTE_PER_CELL};

    //!< Calculate bytes per node (dof-per-stratum)
    IntArrayType<4> node_offset;
    IntScalarType node_counts = 0;
    for (IntScalarType d = 0; d < 4; d++)
    {
        node_offset[d] = node_counts;
        if (bytes_per_stratum[d] > 0)
        {
            node_counts += stratum_ends[d] - stratum_starts[d];
        }
    }

    IntVectorType bytes_per_node(node_counts, 0);
    for (IntScalarType d = 0; d < 4; d++)
    {
        if (bytes_per_stratum[d] == 0)
            continue;

        const IntScalarType& curr_size = stratum_size[d];
        for (IntScalarType i = 0; i < curr_size; i++)
        {
            bytes_per_node[node_offset[d] + i] = bytes_per_stratum[d];
        }
    }

    Int64ScalarType total_bytes =
        accumulate(bytes_per_node.begin( ), bytes_per_node.end( ), Int64ScalarType{0});

    //!<-----------------------------------------------------------------------
    //!< Recover element mesh into arrays for patch partitioning
    //!<-----------------------------------------------------------------------
    IntVectorType elem_cnn(INITIAL_SIZE, -1);
    IntVectorType elem_cnn_offsets(rank_cell_end - rank_cell_start + 1, 0);

    IntScalarType index = 0;
    for (IntScalarType cell_index = rank_cell_start; cell_index < rank_cell_end;
         cell_index++)
    {
        IntScalarType closure_size;
        IntScalarType* closure = NULL;
        DMPlexGetTransitiveClosure(dm, cell_index, true, &closure_size,
                                   &closure);

        for (IntScalarType d = 0; d < 4; d++)
        {
            if (bytes_per_stratum[d] == 0)
                continue;

            for (IntScalarType i = 0; i < 2 * closure_size; i += 2)
            {  //!< Odd index have orientations
                const IntScalarType& curr_point = closure[i];
                IntScalarType depth;
                DMPlexGetPointDepth(dm, curr_point, &depth);
                if (depth == d)
                {
                    if (index >= (IntScalarType)elem_cnn.size( ))
                        elem_cnn.resize(elem_cnn.size( ) + INITIAL_SIZE, -1);

                    //!< Rank local numbering.
                    elem_cnn[index] =
                        node_offset[d] + curr_point - stratum_starts[d];
                    index++;
                }
            }
        }

        elem_cnn_offsets[cell_index - rank_cell_start + 1] = index;
        DMPlexRestoreTransitiveClosure(dm, cell_index, true, &closure_size,
                                       &closure);
    }
    elem_cnn.resize(index);

    //!<-----------------------------------------------------------------------
    //!< Reverse map and adjacency
    //!<-----------------------------------------------------------------------
    IntVectorType node2elem(elem_cnn.size( ), -1);
    IntVectorType node2elem_offsets(node_counts + 1, 0);

    for (IntScalarType n = 0; n < (IntScalarType)elem_cnn.size( ); n++)
        node2elem_offsets[elem_cnn[n] + 1]++;

    for (IntScalarType n = 0; n < node_counts; n++)
        node2elem_offsets[n + 1] += node2elem_offsets[n];

    {
        IntVectorType node2elem_loc(node2elem_offsets.begin( ),
                                    node2elem_offsets.end( ) - 1);
        for (IntScalarType e = 0; e < stratum_size[3]; e++)
        {
            for (IntScalarType n = elem_cnn_offsets[e]; n < elem_cnn_offsets[e + 1];
                 n++)
                node2elem[node2elem_loc[elem_cnn[n]]++] = e;
        }
    }

    //!< Element Adjacency; determined by shared nodes
    //!< AN: These will be small so we can afford push_back.
    IntVectorType elem2elem;
    IntVectorType elem2elem_weights;  //!< Based on shared nodes
    IntVectorType elem2elem_offsets(stratum_size[3] + 1, 0);

    {
        IntVectorType visit_counts(stratum_size[3], 0);
        IntVectorType neighbors;
        for (IntScalarType e = 0; e < stratum_size[3]; e++)
        {
            for (IntScalarType n = elem_cnn_offsets[e]; n < elem_cnn_offsets[e + 1];
                 n++)
            {
                const IntScalarType node = elem_cnn[n];
                for (IntScalarType ee = node2elem_offsets[node];
                     ee < node2elem_offsets[node + 1]; ee++)
                {
                    const IntScalarType curr_elem = node2elem[ee];
                    if (curr_elem == e)
                        continue;

                    if (visit_counts[curr_elem] == 0)
                        neighbors.push_back(curr_elem);
                    visit_counts[curr_elem]++;  //!< += bytes_per_node[node]
                }
            }

            stable_sort(neighbors.begin( ), neighbors.end( ));
            for (IntScalarType ne = 0; ne < (IntScalarType)neighbors.size( ); ne++)
            {
                elem2elem.push_back(neighbors[ne]);
                elem2elem_weights.push_back(visit_counts[neighbors[ne]]);
                visit_counts[neighbors[ne]] = 0;
            }
            neighbors.clear( );
            elem2elem_offsets[e + 1] = elem2elem.size( );
        }
    }

    //!<-----------------------------------------------------------------------
    //!< Compute number of patches
    //!<-----------------------------------------------------------------------

    constexpr RealScalarType WIGGLE_ROOM = 0.75;

    Int64ScalarType byte_upper_bound = cli_data.max_bytes_patch;
    if (byte_upper_bound <= 0 && cli_data.shared_bytes_per_sm > 0)
    {
        byte_upper_bound =
            (cli_data.shared_bytes_per_sm - cli_data.reserve_bytes) /
            cli_data.resident_per_sm;
    }

    IntScalarType npatches;
    IntScalarType npatches_by_bytes = 1;
    IntScalarType npatches_by_elems = 1;
    IntScalarType npatches_by_user = cli_data.npatches;

    if (byte_upper_bound > 0)
    {
        npatches_by_bytes =
            ceil((RealScalarType)total_bytes / (WIGGLE_ROOM * byte_upper_bound));
    }

    if (cli_data.max_elems_patch > 0)
    {
        npatches_by_elems =
            ceil((RealScalarType)stratum_size[3] / cli_data.max_elems_patch);
    }

    const IntScalarType num_slots = cli_data.sm_count * cli_data.resident_per_sm;
    //!< Total patches can be greater than slots available. The patches
    //!< are queued and wait till resources are available.
    const IntScalarType total_gpu_fill = num_slots * cli_data.min_patches_per_sm;

    if (npatches_by_user > 0)
    {
        npatches = npatches_by_user;
    }
    else
    {
        npatches = max({npatches_by_bytes, npatches_by_elems, total_gpu_fill});
        if (num_slots > 0)
            npatches = (npatches + num_slots - 1) / num_slots * num_slots;
    }

    npatches = max(IntScalarType{1}, min(npatches, stratum_size[3]));

    //!<-----------------------------------------------------------------------
    //!< Partition into patches (per GPU block)
    //!<-----------------------------------------------------------------------
    IntVectorType colors(stratum_size[3]);

    if (npatches > 1)
    {
        switch (cli_data.backend)
        {
            case PatchPartitionBackend::METIS:
            {
                const IntScalarType elem2elem_size = elem2elem_offsets[stratum_size[3]];

                PetscSection edge_weights;
                PetscSectionCreate(PETSC_COMM_SELF, &edge_weights);
                PetscSectionSetChart(edge_weights, 0, elem2elem_size);
                for (IntScalarType ee = 0; ee < elem2elem_size; ee++)
                    PetscSectionSetDof(edge_weights, ee, elem2elem_weights[ee]);
                PetscSectionSetUp(edge_weights);

                PetscPartitioner patch_partitioner;
                PetscPartitionerCreate(PETSC_COMM_SELF, &patch_partitioner);
                PetscPartitionerSetType(patch_partitioner,
                                        PETSCPARTITIONERPARMETIS);
                PetscPartitionerSetFromOptions(patch_partitioner);

                PetscSection patch_section;
                IS patch_is;
                PetscSectionCreate(PETSC_COMM_SELF, &patch_section);
                PetscPartitionerPartition(
                    patch_partitioner, npatches, stratum_size[3],
                    elem2elem_offsets.data( ), elem2elem.data( ), NULL,
                    edge_weights, NULL, patch_section, &patch_is);

                const IntScalarType* elem_is = NULL;
                ISGetIndices(patch_is, &elem_is);

                for (IntScalarType p = 0; p < npatches; p++)
                {
                    IntScalarType patch_size, offset;
                    PetscSectionGetDof(patch_section, p, &patch_size);
                    PetscSectionGetOffset(patch_section, p, &offset);

                    for (IntScalarType i = 0; i < patch_size; i++)
                        colors[elem_is[offset + i]] = p;
                }

                ISRestoreIndices(patch_is, &elem_is);

                ISDestroy(&patch_is);
                PetscSectionDestroy(&patch_section);
                PetscSectionDestroy(&edge_weights);
                PetscPartitionerDestroy(&patch_partitioner);

                break;
            }
            case PatchPartitionBackend::JET:
            {
#ifdef HAVE_JET
                LogUtils::PrintError("Jet partitioning not implemented yet.");
#else
                LogUtils::PrintError(
                    "Using Jet paritioner but Jet was not built.");
#endif
            }
        }
    }


    //!<-----------------------------------------------------------------------
    //!< Renumber cells, vertex, edge, and faces based on patches
    //!<-----------------------------------------------------------------------

    //!< Cell order
    IntVectorType elem_per_patch(stratum_size[3]);
    iota(elem_per_patch.begin( ), elem_per_patch.end( ), 0);
    stable_sort(elem_per_patch.begin( ), elem_per_patch.end( ),
                [&](IntScalarType a, IntScalarType b) { return colors[a] < colors[b]; });

    IntVectorType patch_offsets(npatches + 1, 0);
    for (IntScalarType e = 0; e < stratum_size[3]; e++)
        patch_offsets[colors[e] + 1]++;

    for (IntScalarType p = 0; p < npatches; p++)
        patch_offsets[p + 1] += patch_offsets[p];

    elem2patch.resize(stratum_size[3]);
    for (IntScalarType e = 0; e < stratum_size[3]; e++)
        elem2patch[e] = colors[elem_per_patch[e]];  //!< new cell e

    IntVectorType old2new(rank_point_end - rank_point_start, -1);

    //! Renumber cells
    for (IntScalarType e = 0; e < stratum_size[3]; e++)
    {
        old2new[rank_cell_start + elem_per_patch[e] - rank_point_start] =
            rank_cell_start + e;
    }

    //!< Now renumber vertex, edge, and faces
    IntArrayType<3> indices = {rank_vertex_start, rank_edge_start,
                               rank_face_start};
    for (IntScalarType e = 0; e < stratum_size[3]; e++)
    {
        const IntScalarType curr_elem = rank_cell_start + elem_per_patch[e];

        IntScalarType closure_size;
        IntScalarType* closure = NULL;
        DMPlexGetTransitiveClosure(dm, curr_elem, true, &closure_size,
                                   &closure);

        for (IntScalarType k = 0; k < 2 * closure_size; k += 2)
        {
            const IntScalarType curr_point = closure[k];

            //!< Skip the numbered points (first-come-first-serve)
            if (curr_point == rank_cell_start + elem_per_patch[e] ||
                old2new[curr_point - rank_point_start] >= 0)
                continue;

            IntScalarType depth;
            DMPlexGetPointDepth(dm, curr_point, &depth);
            old2new[curr_point - rank_point_start] = indices[depth]++;
        }

        DMPlexRestoreTransitiveClosure(dm, curr_elem, true, &closure_size,
                                       &closure);
    }


    //!<-----------------------------------------------------------------------
    //!< Permute DM (all sections on this DM will follow this order)
    //!<-----------------------------------------------------------------------

    IS old2new_is;
    ISCreateGeneral(PETSC_COMM_SELF, old2new.size( ), old2new.data( ),
                    PETSC_COPY_VALUES, &old2new_is);

    DM dm_perm;
    DMPlexPermute(dm, old2new_is, &dm_perm);
    ISDestroy(&old2new_is);
    DMDestroy(&dm);
    dm = dm_perm;


    //!<-----------------------------------------------------------------------
    //!< Print metrics
    //!<-----------------------------------------------------------------------
    //!< Distinct node bytes per patch (walking cells in patch order)
    Int64ScalarType staged_bytes = 0;
    Int64ScalarType max_patch_bytes = 0;
    IntScalarType min_patch_elems = stratum_size[3];
    IntScalarType max_patch_elems = 0;

    IntVectorType last_patch(node_counts, -1);
    for (IntScalarType p = 0; p < npatches; p++)
    {
        Int64ScalarType patch_bytes = 0;
        for (IntScalarType i = patch_offsets[p]; i < patch_offsets[p + 1]; i++)
        {
            const IntScalarType e = elem_per_patch[i];
            for (IntScalarType n = elem_cnn_offsets[e]; n < elem_cnn_offsets[e + 1];
                 n++)
            {
                const IntScalarType& node = elem_cnn[n];
                if (last_patch[node] != p)
                {
                    last_patch[node] = p;
                    patch_bytes += bytes_per_node[node];
                }
            }
        }

        const IntScalarType& patch_elems = patch_offsets[p + 1] - patch_offsets[p];
        staged_bytes += patch_bytes;
        max_patch_bytes = max(max_patch_bytes, patch_bytes);
        min_patch_elems = min(min_patch_elems, patch_elems);
        max_patch_elems = max(max_patch_elems, patch_elems);
    }

    //!< Connectivity, sum over nodes of (patches touching the node - 1)
    IntScalarType connectivity = 0;
    IntScalarType shared_nodes = 0;

    IntVectorType patch_seen(npatches, -1);  //!< last node that saw patch
    for (IntScalarType node = 0; node < node_counts; node++)
    {
        IntScalarType lambda = 0;
        for (IntScalarType ee = node2elem_offsets[node];
             ee < node2elem_offsets[node + 1]; ee++)
        {
            const IntScalarType p = colors[node2elem[ee]];
            if (patch_seen[p] != node)
            {
                patch_seen[p] = node;
                lambda++;
            }
        }

        if (lambda > 1)
        {
            connectivity += lambda - 1;
            shared_nodes++;
        }
    }

    //!< Reduce to rank 0
    IntScalarType g_npatches_sum, g_npatches_min, g_npatches_max;
    IntScalarType g_min_patch_elems, g_max_patch_elems;
    IntScalarType g_connectivity, g_shared_nodes, g_nodes;
    Int64ScalarType g_total_bytes, g_staged_bytes, g_max_patch_bytes;

    MPI_Reduce(&npatches, &g_npatches_sum, 1, MPIU_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&npatches, &g_npatches_min, 1, MPIU_INT, MPI_MIN, 0, comm);
    MPI_Reduce(&npatches, &g_npatches_max, 1, MPIU_INT, MPI_MAX, 0, comm);
    MPI_Reduce(&min_patch_elems, &g_min_patch_elems, 1, MPIU_INT, MPI_MIN, 0,
               comm);
    MPI_Reduce(&max_patch_elems, &g_max_patch_elems, 1, MPIU_INT, MPI_MAX, 0,
               comm);
    MPI_Reduce(&connectivity, &g_connectivity, 1, MPIU_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&shared_nodes, &g_shared_nodes, 1, MPIU_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&node_counts, &g_nodes, 1, MPIU_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&total_bytes, &g_total_bytes, 1, MPI_INT64_T, MPI_SUM, 0, comm);
    MPI_Reduce(&staged_bytes, &g_staged_bytes, 1, MPI_INT64_T, MPI_SUM, 0,
               comm);
    MPI_Reduce(&max_patch_bytes, &g_max_patch_bytes, 1, MPI_INT64_T, MPI_MAX, 0,
               comm);

    //!< Only rank 0 values are meaningful (and printed)
    LogUtils::Print("Patch metrics (all ranks):");
    LogUtils::Print("  patches (total / min / max per rank) : %d / %d / %d",
                    g_npatches_sum, g_npatches_min, g_npatches_max);
    LogUtils::Print("  elements per patch (min / max)       : %d / %d",
                    g_min_patch_elems, g_max_patch_elems);
    LogUtils::Print("  max bytes per patch                  : %lld",
                    (long long)g_max_patch_bytes);
    LogUtils::Print("  total / staged bytes                 : %lld / %lld",
                    (long long)g_total_bytes, (long long)g_staged_bytes);
    LogUtils::Print(
        "  replication (staged / total)         : %.3f",
        g_total_bytes > 0 ? (double)g_staged_bytes / g_total_bytes : 0.0);
    LogUtils::Print("  connectivity (sum lambda - 1)        : %d",
                    g_connectivity);
    LogUtils::Print("  shared node fraction                 : %.3f",
                    g_nodes > 0 ? (double)g_shared_nodes / g_nodes : 0.0);
}

//-----------------------------------------------------------------------------

void SpaceDiscretization::Destroy( )
{
    DMDestroy(&dm);
}

//-----------------------------------------------------------------------------

void SpaceDiscretization::OutputMesh( )
{
    Int32ScalarType rank;
    MPI_Comm_rank(comm, &rank);

    //!<-----------------------------------------------------------------------
    //!< Step 1: Create cell-wise color vector
    //!<-----------------------------------------------------------------------
    //!< Create a (ndof, offset) tuple for the colors
    PetscSection section;
    PetscSectionCreate(comm, &section);
    //!< Without fields we get "color(null)" in ParaView; not sure why.
    PetscSectionSetNumFields(section, 2);
    PetscSectionSetFieldName(section, 0, "rank");
    PetscSectionSetFieldName(section, 1, "patch");
    PetscSectionSetChart(section, rank_point_start, rank_point_end);
    for (IntScalarType cell_index = rank_cell_start; cell_index < rank_cell_end;
         cell_index++)
    {
        PetscSectionSetDof(section, cell_index, 2);
        PetscSectionSetFieldDof(section, cell_index, 0, 1);
        PetscSectionSetFieldDof(section, cell_index, 1, 1);
    }
    PetscSectionSetUp(section);

    //!< Create a duplicate DM; used only for mesh visualization
    DM dm_colors;
    DMClone(dm, &dm_colors);
    DMSetLocalSection(dm_colors, section);

    //!<-----------------------------------------------------------------------
    //!< Step 2: Create local vector
    //!<-----------------------------------------------------------------------
    RealScalarType* colors_data;
    Vec colors;

    DMCreateLocalVector(dm_colors, &colors);
    PetscObjectSetName((PetscObject)colors, "");

    //!<-----------------------------------------------------------------------
    //!< Step 3: Fill local vector with ranks
    //!<-----------------------------------------------------------------------
    //!< Shift local patch ids so every patch is unique across ranks
    IntScalarType rank_npatches = 0;
    for (IntScalarType e = 0; e < (IntScalarType)elem2patch.size( ); e++)
        rank_npatches = max(rank_npatches, elem2patch[e] + 1);

    IntScalarType patch_shift = 0;
    MPI_Exscan(&rank_npatches, &patch_shift, 1, MPIU_INT, MPI_SUM, comm);
    if (rank == 0)
        patch_shift = 0;  //!< Exscan leaves rank 0 undefined

    VecGetArray(colors, &colors_data);

    for (IntScalarType cell_index = rank_cell_start; cell_index < rank_cell_end;
         cell_index++)
    {
        for (IntScalarType field_index = 0; field_index < 2; field_index++)
        {
            IntScalarType num_dofs, offset;
            PetscSectionGetFieldDof(section, cell_index, field_index,
                                    &num_dofs);
            PetscSectionGetFieldOffset(section, cell_index, field_index,
                                       &offset);

            for (IntScalarType dof = 0; dof < num_dofs; dof++)
                colors_data[offset + dof] =
                    (field_index == 0)
                        ? rank
                        : patch_shift +
                              elem2patch[cell_index - rank_cell_start];
        }
    }

    VecRestoreArray(colors, &colors_data);

    //!<-----------------------------------------------------------------------
    //!< Step 4: Create viewer and visualize mesh
    //!<-----------------------------------------------------------------------
    PetscViewer viewer;
    PetscViewerVTKOpen(comm, cli_data.out_file.c_str( ), FILE_MODE_WRITE,
                       &viewer);
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
