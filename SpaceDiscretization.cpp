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
#include <cstdlib>
#include <limits>
#include <numeric>

using std::abs;
using std::accumulate;
using std::ceil;
using std::floor;
using std::iota;
using std::log2;
using std::max;
using std::min;
using std::pow;
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

    Int64VectorType bytes_per_node(node_counts, 0);
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

    Int64ScalarType total_bytes = accumulate(
        bytes_per_node.begin( ), bytes_per_node.end( ), Int64ScalarType{0});

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
            for (IntScalarType n = elem_cnn_offsets[e];
                 n < elem_cnn_offsets[e + 1]; n++)
                node2elem[node2elem_loc[elem_cnn[n]]++] = e;
        }
    }

    //!< Element Adjacency; determined by shared nodes
    //!< AN: These will be small so we can afford push_back.
    RealVectorType elem_weights(stratum_size[3]);  //!< weights on each element
    IntVectorType elem2elem;
    IntVectorType elem2elem_weights;  //!< elem-elem edge weights
    IntVectorType elem2elem_offsets(stratum_size[3] + 1, 0);

    {
        IntIntScalarMap neighbor_weights;
        for (IntScalarType e = 0; e < stratum_size[3]; e++)
        {
            RealScalarType wgt = 0;
            for (IntScalarType n = elem_cnn_offsets[e];
                 n < elem_cnn_offsets[e + 1]; n++)
            {
                const IntScalarType node = elem_cnn[n];

                wgt += (RealScalarType)bytes_per_node[node] /
                       (node2elem_offsets[node + 1] - node2elem_offsets[node]);

                for (IntScalarType ee = node2elem_offsets[node];
                     ee < node2elem_offsets[node + 1]; ee++)
                {
                    const IntScalarType curr_elem = node2elem[ee];
                    if (curr_elem == e)
                        continue;

                    neighbor_weights[curr_elem]++;  //!< += bytes_per_node[node]
                }
            }

            elem_weights[e] = wgt;
            for (IISMCIter it = neighbor_weights.begin( );
                 it != neighbor_weights.end( ); ++it)
            {
                elem2elem.push_back(it->first);
                elem2elem_weights.push_back(it->second);
            }
            neighbor_weights.clear( );
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
            cli_data.shared_bytes_per_sm / cli_data.resident_per_sm -
            cli_data.reserve_bytes;  //!< reserve is per patch
    }

    IntScalarType npatches;
    IntScalarType npatches_by_bytes = 1;
    IntScalarType npatches_by_elems = 1;
    IntScalarType npatches_by_user = cli_data.npatches;

    if (byte_upper_bound > 0)
    {
        npatches_by_bytes = ceil((RealScalarType)total_bytes /
                                 (WIGGLE_ROOM * byte_upper_bound));
    }

    if (cli_data.max_elems_patch > 0)
    {
        npatches_by_elems =
            ceil((RealScalarType)stratum_size[3] / cli_data.max_elems_patch);
    }

    const IntScalarType num_slots =
        cli_data.sm_count * cli_data.resident_per_sm;
    const IntScalarType total_gpu_fill =
        cli_data.sm_count * cli_data.min_patches_per_sm;

    if (npatches_by_user > 0)
    {
        npatches = npatches_by_user;
    }
    else
    {
        npatches = max({npatches_by_bytes, npatches_by_elems, total_gpu_fill});
        if (num_slots > 0)
            npatches = ceil((RealScalarType)npatches / num_slots) * num_slots;
    }

    npatches = max(IntScalarType{1}, min(npatches, stratum_size[3]));

    //!< Do nothing for a single patch.
    if (npatches == IntScalarType{1})
    {
        return;
    }

    //!<-----------------------------------------------------------------------
    //!< Partition into patches (per GPU block)
    //!<-----------------------------------------------------------------------
    IntVectorType colors(stratum_size[3]);

    //!< node_marker[n] == curr_pass: node n already counted in this pass
    IntVectorType node_marker(node_counts, -1);
    IntScalarType curr_pass = -1;

    switch (cli_data.backend)
    {
        case PatchPartitionBackend::METIS:
        {
            const IntScalarType elem2elem_size =
                elem2elem_offsets[stratum_size[3]];

            PetscSection edge_weights;
            PetscSectionCreate(PETSC_COMM_SELF, &edge_weights);
            PetscSectionSetChart(edge_weights, 0, elem2elem_size);
            for (IntScalarType ee = 0; ee < elem2elem_size; ee++)
            {
                PetscSectionSetDof(edge_weights, ee, elem2elem_weights[ee]);
            }
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

            PetscSectionDestroy(&edge_weights);
            PetscPartitionerDestroy(&patch_partitioner);

            //!< Store patch based colors
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

            //!< Finalize and Destroy temporary objects
            ISRestoreIndices(patch_is, &elem_is);
            ISDestroy(&patch_is);
            PetscSectionDestroy(&patch_section);

            break;
        }
        case PatchPartitionBackend::JET:
        {
#ifdef HAVE_JET
            LogUtils::PrintError("Jet partitioning not implemented yet.");
#else
            LogUtils::PrintError("Using Jet paritioner but Jet was not built.");
#endif
            break;
        }
        case PatchPartitionBackend::CCHP:
        {
            if (byte_upper_bound <= 0)
            {
                LogUtils::PrintError(
                    "CCHP needs --bytes-per-patch or "
                    "--shared-bytes-per-sm.");
                break;
            }
            if (npatches_by_user > 0)
            {
                LogUtils::Print(
                    "CCHP sets the patch count itself; --npatches is "
                    "ignored.");
            }

            //!< Fewest patches per rank, rounded up to whole GPU rounds
            IntScalarType min_patches = max(npatches_by_elems, total_gpu_fill);
            if (num_slots > 0)
            {
                min_patches =
                    ceil((RealScalarType)min_patches / num_slots) * num_slots;
            }
            min_patches =
                max(IntScalarType{1}, min(min_patches, stratum_size[3]));

            const IntScalarType max_elems = (cli_data.max_elems_patch > 0)
                                                ? cli_data.max_elems_patch
                                                : stratum_size[3];
            constexpr RealScalarType IMBALANCE = 0.05;

            //!< Initialize
            PetscPartitioner bisector;
            PetscPartitionerCreate(PETSC_COMM_SELF, &bisector);
            PetscPartitionerSetType(bisector, PETSCPARTITIONERPARMETIS);
            PetscPartitionerSetFromOptions(bisector);

            IntVector2DType patches_dfs;  //!< All patches in DFS order
            Int64VectorType patch_bytes;
            IntVectorType patch2rank(stratum_size[3], -1);
            IntVectorType patch_side(stratum_size[3], -1);

            //!< node_count[s][n]: elements on side s touching node n
            IntVector2DType node_count(2, IntVectorType(node_counts, 0));

            //!< Temporary stacks: (part, patches this part should become)
            IntStack2DType patch_stack;
            IntStackType split_target;

            //!< Push all elements in the stack
            patch_stack.push(IntVectorType(stratum_size[3]));
            iota(patch_stack.top( ).begin( ), patch_stack.top( ).end( ), 0);
            split_target.push(min_patches);

            while (!patch_stack.empty( ))
            {
                const IntScalarType target = split_target.top( );
                IntVectorType& curr_patch = patch_stack.top( );
                split_target.pop( );

                const IntScalarType curr_patch_size = curr_patch.size( );

                //!<-----------------------------------------------------------
                //!< Step 1a: Calculate exact staged bytes for this partition
                //!<-----------------------------------------------------------
                Int64ScalarType curr_patch_bytes = 0;
                curr_pass++;
                for (IntScalarType e = 0; e < curr_patch_size; e++)
                {
                    const IntScalarType curr_elem = curr_patch[e];
                    for (IntScalarType n = elem_cnn_offsets[curr_elem];
                         n < elem_cnn_offsets[curr_elem + 1]; n++)
                    {
                        const IntScalarType node = elem_cnn[n];
                        if (node_marker[node] != curr_pass)
                        {
                            node_marker[node] = curr_pass;
                            curr_patch_bytes += bytes_per_node[node];
                        }
                    }
                }

                //!< Patches needed: by bytes, by elements, or its target
                const IntScalarType capacity_ratio = max(
                    {(IntScalarType)ceil((RealScalarType)curr_patch_bytes /
                                         byte_upper_bound),
                     (IntScalarType)ceil((RealScalarType)curr_patch_size /
                                         max_elems),
                     target});

                //!<-----------------------------------------------------------
                //!< Step 1b: This patch can fit
                //!<-----------------------------------------------------------
                if (capacity_ratio <= 1 || curr_patch_size == 1)
                {
                    //!< AN: This might be expensive. Will improve later.
                    patches_dfs.push_back(curr_patch);
                    patch_bytes.push_back(curr_patch_bytes);
                    patch_stack.pop( );
                    continue;
                }

                //!<-----------------------------------------------------------
                //!< Step 2a: Split partition into two parts
                //!<-----------------------------------------------------------
                const IntScalarType left_split_ratio = capacity_ratio / 2;
                const IntScalarType right_split_ratio =
                    capacity_ratio - left_split_ratio;

                for (IntScalarType e = 0; e < curr_patch_size; e++)
                {
                    patch2rank[curr_patch[e]] = e;
                }

                IntVectorType patch_adj;
                IntVectorType patch_adj_weights;
                IntVectorType patch_adj_offsets(curr_patch_size + 1, 0);

                for (IntScalarType e = 0; e < curr_patch_size; e++)
                {
                    const IntScalarType curr_elem = curr_patch[e];
                    for (IntScalarType k = elem2elem_offsets[curr_elem];
                         k < elem2elem_offsets[curr_elem + 1]; k++)
                    {
                        const IntScalarType f = patch2rank[elem2elem[k]];
                        if (f >= 0)
                        {
                            patch_adj.push_back(f);
                            patch_adj_weights.push_back(elem2elem_weights[k]);
                        }
                    }
                    patch_adj_offsets[e + 1] = patch_adj.size( );
                }

                for (IntScalarType e = 0; e < curr_patch_size; e++)
                {
                    patch2rank[curr_patch[e]] = -1;
                }

                //!<-----------------------------------------------------------
                //!< Step 2b: Split using METIS
                //!<-----------------------------------------------------------
                //!< No vertex weights: METIS balances element count
                PetscSection edge_section, target_section, split_section;

                PetscSectionCreate(PETSC_COMM_SELF, &edge_section);
                PetscSectionSetChart(edge_section, 0, patch_adj.size( ));
                for (IntScalarType k = 0; k < (IntScalarType)patch_adj.size( );
                     k++)
                {
                    PetscSectionSetDof(edge_section, k, patch_adj_weights[k]);
                }
                PetscSectionSetUp(edge_section);

                PetscSectionCreate(PETSC_COMM_SELF, &target_section);
                PetscSectionSetChart(target_section, 0, 2);
                PetscSectionSetDof(target_section, 0, left_split_ratio);
                PetscSectionSetDof(target_section, 1, right_split_ratio);
                PetscSectionSetUp(target_section);

                IS split_is;
                PetscSectionCreate(PETSC_COMM_SELF, &split_section);
                PetscPartitionerPartition(
                    bisector, 2, curr_patch_size, patch_adj_offsets.data( ),
                    patch_adj.data( ), NULL, edge_section,
                    target_section, split_section, &split_is);

                IntScalarType num_left;
                const IntScalarType* split_ids = NULL;
                PetscSectionGetDof(split_section, 0, &num_left);
                ISGetIndices(split_is, &split_ids);
                for (IntScalarType e = 0; e < curr_patch_size; e++)
                {
                    patch_side[curr_patch[split_ids[e]]] =
                        (e < num_left) ? 0 : 1;
                }
                ISRestoreIndices(split_is, &split_ids);

                ISDestroy(&split_is);
                PetscSectionDestroy(&split_section);
                PetscSectionDestroy(&target_section);
                PetscSectionDestroy(&edge_section);

                if (num_left == 0 || num_left == curr_patch_size)
                {
                    const IntScalarType cut =
                        max(IntScalarType{1}, curr_patch_size *
                                                  left_split_ratio /
                                                  capacity_ratio);
                    for (IntScalarType e = 0; e < curr_patch_size; e++)
                    {
                        patch_side[curr_patch[e]] = (e < cut) ? 0 : 1;
                    }
                }

                //!<-----------------------------------------------------------
                //!< Step 3: Refinement and corrections
                //!<-----------------------------------------------------------
                IntArrayType<2> side_sizes = {0, 0};

                for (IntScalarType e = 0; e < curr_patch_size; e++)
                {
                    const IntScalarType curr_elem = curr_patch[e];
                    const IntScalarType s = patch_side[curr_elem];
                    side_sizes[s]++;
                    for (IntScalarType n = elem_cnn_offsets[curr_elem];
                         n < elem_cnn_offsets[curr_elem + 1]; n++)
                    {
                        node_count[s][elem_cnn[n]]++;
                    }
                }

                //!< Target element count of each side, and the most it holds
                const RealArrayType<2> target_size = {
                    (RealScalarType)curr_patch_size * left_split_ratio /
                        capacity_ratio,
                    (RealScalarType)curr_patch_size * right_split_ratio /
                        capacity_ratio};

                //!< Spread IMBALANCE over the cuts below; avoids compounding
                const IntScalarType levels =
                    ceil(log2((RealScalarType)capacity_ratio));
                const RealScalarType slack =
                    pow(1.0 + IMBALANCE, 1.0 / levels) - 1.0;

                const IntArrayType<2> size_limit = {
                    (IntScalarType)ceil((1 + slack) * target_size[0]),
                    (IntScalarType)ceil((1 + slack) * target_size[1])};

                //!< Balance first: while a side is over its limit, move out
                //!< the element with the highest (possibly negative) gain.
                //!< Limits sum to >= curr_patch_size, so at most one side is
                //!< over, and the other side cannot go over its own limit.
                for (IntScalarType from = 0; from < 2; from++)
                {
                    const IntScalarType to = 1 - from;
                    while (side_sizes[from] > size_limit[from])
                    {
                        IntScalarType best_elem = -1;
                        Int64ScalarType best_gain =
                            std::numeric_limits<Int64ScalarType>::lowest( );

                        //!< Fixed order; first best wins (deterministic)
                        for (IntScalarType e = 0; e < curr_patch_size; e++)
                        {
                            const IntScalarType curr_elem = curr_patch[e];
                            if (patch_side[curr_elem] != from)
                                continue;

                            Int64ScalarType removed = 0;
                            Int64ScalarType added = 0;
                            for (IntScalarType n = elem_cnn_offsets[curr_elem];
                                 n < elem_cnn_offsets[curr_elem + 1]; n++)
                            {
                                const IntScalarType node = elem_cnn[n];
                                if (node_count[from][node] == 1)
                                    removed += bytes_per_node[node];
                                if (node_count[to][node] == 0)
                                    added += bytes_per_node[node];
                            }

                            if (removed - added > best_gain)
                            {
                                best_gain = removed - added;
                                best_elem = curr_elem;
                            }
                        }

                        for (IntScalarType n = elem_cnn_offsets[best_elem];
                             n < elem_cnn_offsets[best_elem + 1]; n++)
                        {
                            const IntScalarType node = elem_cnn[n];
                            node_count[from][node]--;
                            node_count[to][node]++;
                        }

                        side_sizes[from]--;
                        side_sizes[to]++;
                        patch_side[best_elem] = to;
                    }
                }

                bool moved = true;
                while (moved)
                {
                    moved = false;
                    for (IntScalarType e = 0; e < curr_patch_size; e++)
                    {
                        const IntScalarType curr_elem = curr_patch[e];
                        const IntScalarType from = patch_side[curr_elem];
                        const IntScalarType to = 1 - from;

                        Int64ScalarType removed = 0;
                        Int64ScalarType added = 0;
                        for (IntScalarType n = elem_cnn_offsets[curr_elem];
                             n < elem_cnn_offsets[curr_elem + 1]; n++)
                        {
                            const IntScalarType node = elem_cnn[n];
                            if (node_count[from][node] == 1)
                                removed += bytes_per_node[node];
                            if (node_count[to][node] == 0)
                                added += bytes_per_node[node];
                        }

                        const Int64ScalarType gain = removed - added;

                        //!< Distance from the target split, before and after
                        const RealScalarType before =
                            abs(side_sizes[from] - target_size[from]) +
                            abs(side_sizes[to] - target_size[to]);
                        const RealScalarType after =
                            abs(side_sizes[from] - 1 - target_size[from]) +
                            abs(side_sizes[to] + 1 - target_size[to]);

                        if (!(gain > 0 || (gain == 0 && after < before)))
                            continue;
                        if (side_sizes[to] + 1 > size_limit[to] ||
                            side_sizes[from] == 1)
                            continue;

                        for (IntScalarType n = elem_cnn_offsets[curr_elem];
                             n < elem_cnn_offsets[curr_elem + 1]; n++)
                        {
                            const IntScalarType node = elem_cnn[n];
                            node_count[from][node]--;
                            node_count[to][node]++;
                        }

                        side_sizes[from]--;
                        side_sizes[to]++;
                        patch_side[curr_elem] = to;
                        moved = true;
                    }
                }

                //!<-----------------------------------------------------------
                //!< Step 4: Push right, then left, each with its share
                //!<-----------------------------------------------------------
                IntVectorType left_patch;
                IntVectorType right_patch;
                for (IntScalarType e = 0; e < curr_patch_size; e++)
                {
                    //!< Reset counts touched by this part for the next split
                    const IntScalarType curr_elem = curr_patch[e];
                    for (IntScalarType n = elem_cnn_offsets[curr_elem];
                         n < elem_cnn_offsets[curr_elem + 1]; n++)
                    {
                        node_count[0][elem_cnn[n]] = 0;
                        node_count[1][elem_cnn[n]] = 0;
                    }

                    if (patch_side[curr_patch[e]] == 0)
                        left_patch.push_back(curr_patch[e]);
                    else
                        right_patch.push_back(curr_patch[e]);
                }

                patch_stack.pop( );
                patch_stack.push(right_patch);
                split_target.push(right_split_ratio);
                patch_stack.push(left_patch);
                split_target.push(left_split_ratio);
            }

            npatches = patches_dfs.size( );
            for (IntScalarType p = 0; p < npatches; p++)
            {
                for (IntScalarType i = 0;
                     i < (IntScalarType)patches_dfs[p].size( ); i++)
                    colors[patches_dfs[p][i]] = p;
            }

            PetscPartitionerDestroy(&bisector);
            break;
        }
    }

    //!<-----------------------------------------------------------------------
    //!< Renumber cells, vertex, edge, and faces based on patches
    //!<-----------------------------------------------------------------------

    //!< Cell order
    IntVectorType elem_per_patch(stratum_size[3]);
    iota(elem_per_patch.begin( ), elem_per_patch.end( ), 0);
    stable_sort(elem_per_patch.begin( ), elem_per_patch.end( ),
                [&](IntScalarType a, IntScalarType b)
                { return colors[a] < colors[b]; });

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

    //!< DMPlexPermute does not carry the point SF; rebuild it
    PetscSF sf_old;
    DMGetPointSF(dm, &sf_old);

    IntScalarType nroots, nleaves;
    const IntScalarType* ilocal = NULL;
    const PetscSFNode* iremote = NULL;
    PetscSFGetGraph(sf_old, &nroots, &nleaves, &ilocal, &iremote);

    //!< nroots < 0: graph never set (serial run, no DMPlexDistribute)
    if (nroots >= 0)
    {
        //!< New number of each leaf's owner point, as numbered on the owner
        IntVectorType remote_new(old2new.size( ), -1);
        PetscSFBcastBegin(sf_old, MPIU_INT, old2new.data( ), remote_new.data( ),
                          MPI_REPLACE);
        PetscSFBcastEnd(sf_old, MPIU_INT, old2new.data( ), remote_new.data( ),
                        MPI_REPLACE);

        IntScalarType* new_local;
        PetscSFNode* new_remote;
        PetscMalloc1(nleaves, &new_local);
        PetscMalloc1(nleaves, &new_remote);
        for (IntScalarType l = 0; l < nleaves; l++)
        {
            const IntScalarType p = ilocal ? ilocal[l] : l;
            new_local[l] = old2new[p - rank_point_start];
            new_remote[l].rank = iremote[l].rank;
            new_remote[l].index = remote_new[p - rank_point_start];
        }

        PetscSF sf_new;
        PetscSFCreate(comm, &sf_new);
        PetscSFSetGraph(sf_new, nroots, nleaves, new_local, PETSC_OWN_POINTER,
                        new_remote, PETSC_OWN_POINTER);
        DMSetPointSF(dm_perm, sf_new);

        //!< Coordinate DM keeps its own SF (not propagated by DMSetPointSF)
        DM cdm_perm;
        DMGetCoordinateDM(dm_perm, &cdm_perm);
        DMSetPointSF(cdm_perm, sf_new);
        PetscSFDestroy(&sf_new);
    }

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

    for (IntScalarType p = 0; p < npatches; p++)
    {
        Int64ScalarType patch_bytes = 0;
        curr_pass++;
        for (IntScalarType i = patch_offsets[p]; i < patch_offsets[p + 1]; i++)
        {
            const IntScalarType e = elem_per_patch[i];
            for (IntScalarType n = elem_cnn_offsets[e];
                 n < elem_cnn_offsets[e + 1]; n++)
            {
                const IntScalarType& node = elem_cnn[n];
                if (node_marker[node] != curr_pass)
                {
                    node_marker[node] = curr_pass;
                    patch_bytes += bytes_per_node[node];
                }
            }
        }

        const IntScalarType& patch_elems =
            patch_offsets[p + 1] - patch_offsets[p];
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
        const IntScalarType patch_id =
            patch_shift + elem2patch[cell_index - rank_cell_start];

        for (IntScalarType field_index = 0; field_index < 2; field_index++)
        {
            IntScalarType num_dofs, offset;
            PetscSectionGetFieldDof(section, cell_index, field_index,
                                    &num_dofs);
            PetscSectionGetFieldOffset(section, cell_index, field_index,
                                       &offset);

            const RealScalarType color = (field_index == 0) ? rank : patch_id;
            for (IntScalarType dof = 0; dof < num_dofs; dof++)
                colors_data[offset + dof] = color;
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
