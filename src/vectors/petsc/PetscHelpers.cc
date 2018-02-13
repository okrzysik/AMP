// This file contains so definitions and wrapper functions for PETSc
#ifndef PETSC_HELPERS
#define PETSC_HELPERS

#include "AMP/vectors/petsc/PetscHelpers.h"

#include "petsc.h"
#include "petscvec.h"


namespace PETSC {


/********************************************************
 * Helper functions                                      *
 ********************************************************/
PetscErrorCode vecDestroy( Vec *v )
{
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    return VecDestroy( *v );
#elif PETSC_VERSION_GE( 3, 2, 0 )
    return VecDestroy( v );
#else
#error Not programmed for this version yet
#endif
}
PetscErrorCode randomDestroy( PetscRandom *random )
{
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    return PetscRandomDestroy( *random );
#elif PETSC_VERSION_GE( 3, 2, 0 )
    return PetscRandomDestroy( random );
#else
#error Not programmed for this version of petsc
#endif
}
PetscErrorCode matDestroy( Mat *mat )
{
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    return MatDestroy( *mat );
#elif PETSC_VERSION_GE( 3, 2, 0 )
    return MatDestroy( mat );
#else
#error Not programmed for this version yet
#endif
}


} // namespace PETSC

#endif
