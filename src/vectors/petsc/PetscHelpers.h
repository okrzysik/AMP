// This file contains so definitions and wrapper functions for PETSc
#ifndef PETSC_HELPERS
#define PETSC_HELPERS

#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/enable_shared_from_this.h"
#include "AMP/vectors/data/DataChangeListener.h"

#include <memory>


// Forward declare a few types for PETSc
typedef int PetscErrorCode;
typedef struct _p_Vec *Vec;
typedef struct _p_Mat *Mat;
typedef struct _p_PetscRandom *PetscRandom;


namespace AMP::LinearAlgebra {
class Vector;
class Matrix;
} // namespace AMP::LinearAlgebra

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


/********************************************************
 * Get the AMP vector from the PETSc Vec or Mat          *
 ********************************************************/
std::shared_ptr<AMP::LinearAlgebra::Vector> getAMP( Vec t );
std::shared_ptr<AMP::LinearAlgebra::Matrix> getAMP( Mat t );


/********************************************************
 * Get a PETSc Vec or Mat from an AMP vector/matrix      *
 ********************************************************/
Vec getVec( std::shared_ptr<AMP::LinearAlgebra::Vector> v );
Mat getMat( std::shared_ptr<AMP::LinearAlgebra::Matrix> m );


} // namespace PETSC

#endif
