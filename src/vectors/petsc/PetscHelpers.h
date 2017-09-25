// This file contains so definitions and wrapper functions for PETSc
#ifndef PETSC_HELPERS
#define PETSC_HELPERS

#include "petscmat.h"
#include "petscvec.h"


#ifndef PETSC_VERSION_LT

#define PETSC_VERSION_LT( MAJOR, MINOR, SUBMINOR )                              \
    ( PETSC_VERSION_RELEASE == 1 && ( PETSC_VERSION_MAJOR < ( MAJOR ) ||        \
                                      ( PETSC_VERSION_MAJOR == ( MAJOR ) &&     \
                                        ( PETSC_VERSION_MINOR < ( MINOR ) ||    \
                                          ( PETSC_VERSION_MINOR == ( MINOR ) && \
                                            ( PETSC_VERSION_SUBMINOR < ( SUBMINOR ) ) ) ) ) ) )

#define PETSC_VERSION_LE( MAJOR, MINOR, SUBMINOR ) \
    ( PETSC_VERSION_LT( MAJOR, MINOR, SUBMINOR ) || PETSC_VERSION_( MAJOR, MINOR, SUBMINOR ) )

#define PETSC_VERSION_GT( MAJOR, MINOR, SUBMINOR ) \
    ( 0 == PETSC_VERSION_LE( MAJOR, MINOR, SUBMINOR ) )

#define PETSC_VERSION_GE( MAJOR, MINOR, SUBMINOR ) \
    ( 0 == PETSC_VERSION_LT( MAJOR, MINOR, SUBMINOR ) )

#endif


namespace PETSC {


/********************************************************
 * Helper functions                                      *
 ********************************************************/
inline PetscErrorCode vecDestroy( Vec *v )
{
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    return VecDestroy( *v );
#elif PETSC_VERSION_GE( 3, 2, 0 )
    return VecDestroy( v );
#else
#error Not programmed for this version yet
#endif
}
inline PetscErrorCode randomDestroy( PetscRandom *random )
{
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    return PetscRandomDestroy( *random );
#elif PETSC_VERSION_GE( 3, 2, 0 )
    return PetscRandomDestroy( random );
#else
#error Not programmed for this version of petsc
#endif
}
inline PetscErrorCode matDestroy( Mat *mat )
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
