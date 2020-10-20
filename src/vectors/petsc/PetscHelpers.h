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


namespace AMP::LinearAlgebra {
class Vector;
class Matrix;
class ManagedPetscVector;
class ManagedPetscMatrix;
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
 * Wrapper class for an AMP vector for PETSc Vec         *
 ********************************************************/
class PetscVectorWrapper
{
public:
    PetscVectorWrapper()                             = delete;
    PetscVectorWrapper( const PetscVectorWrapper & ) = delete;
    PetscVectorWrapper( AMP::LinearAlgebra::ManagedPetscVector *vec );
    ~PetscVectorWrapper();
    bool petscHoldsView() const;
    inline Vec &getVec() { return d_petscVec; }
    inline auto getAMP() { return d_vec; }
    bool check() const;
    Vec duplicate();
    void destroy();

protected:
    Vec d_petscVec;
    bool d_madeWithPetscDuplicate;
    uint32_t hash;
    AMP::LinearAlgebra::ManagedPetscVector *d_vec;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_vec2;
};


} // namespace PETSC

#endif
