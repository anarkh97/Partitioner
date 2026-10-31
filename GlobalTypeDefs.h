#ifndef __GLOBAL_TYPE_DEFS_H__
#define __GLOBAL_TYPE_DEFS_H__

#include <petscsys.h>
#include <array>
#include <map>
#include <set>
#include <string>
#include <vector>

/**************************************
* Global Type definitions
* Note: Most of Kokkos/METIS types and
* other host vs. device types (like,
* FP16, FP32, etc.) are already handled
* by PETSc (see their petscsystypes.h)
**************************************/

// AN: We follow,
// <Ordinal><DataStructure>Type
// naming scheme.

//!< Scalars
using StringScalarType = std::string;
using RealScalarType = PetscReal;
using IntScalarType = PetscInt;
using Int32ScalarType = std::int32_t;
using Int64ScalarType = std::int64_t;

//!< Vectors (dynamic memory allocation)
using StringVectorType = std::vector<StringScalarType>;
using RealVectorType = std::vector<RealScalarType>;
using IntVectorType = std::vector<IntScalarType>;
using Int32VectorType = std::vector<Int32ScalarType>;
using Int64VectorType = std::vector<Int64ScalarType>;

//!< Vectors (static memory allocation)
template<IntScalarType N>
using StringArrayType = std::array<StringScalarType, N>;
template<IntScalarType N>
using RealArrayType = std::array<RealScalarType, N>;
template<IntScalarType N>
using IntArrayType = std::array<IntScalarType, N>;
template<IntScalarType N>
using Int32ArrayType = std::array<Int32ScalarType, N>;
template<IntScalarType N>
using Int64ArrayType = std::array<Int64ScalarType, N>;

//!< Sets
using StringSetType = std::set<StringScalarType>;
using RealSetType = std::set<RealScalarType>;
using IntSetType = std::set<IntScalarType>;
using Int32SetType = std::set<Int32ScalarType>;
using Int64SetType = std::set<Int64ScalarType>;

//!< Maps

#endif
