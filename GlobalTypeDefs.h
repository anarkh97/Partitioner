#ifndef __GLOBAL_TYPE_DEFS_H__
#define __GLOBAL_TYPE_DEFS_H__

#include <petscsys.h>

/**************************************
* Global Type definitions
* Note: Most of Kokkos/METIS types and
* other host vs. device types (like,
* FP16, FP32, etc.) are already handled
* by PETSc (see their petscsystypes.h)
**************************************/

typedef PetscReal real_t;
typedef PetscInt int_t;
typedef std::int64_t byte_t;

#endif
