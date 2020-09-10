// This file contains so definitions and wrapper functions for PETSc
#ifndef PETSC_HELPERS
#define PETSC_HELPERS

#include "AMP/utils/AMP_MPI.h"

#include <memory>


// Forward declare a few types for PETSc
typedef int PetscErrorCode;
typedef struct _p_Vec *Vec;
typedef struct _p_Mat *Mat;
typedef struct _p_PetscRandom *PetscRandom;


namespace PETSC {


/********************************************************
 * Destructors                                           *
 ********************************************************/
PetscErrorCode vecDestroy( Vec *v );
PetscErrorCode randomDestroy( PetscRandom *random );
PetscErrorCode matDestroy( Mat *mat );


/********************************************************
 * Create random number generator                        *
 ********************************************************/
std::shared_ptr<PetscRandom> genPetscRandom( const AMP::AMP_MPI &comm );


/********************************************************
 * Reset petsc vector operations                          *
 ********************************************************/
void reset_vec_ops( Vec t );


} // namespace PETSC

#endif
