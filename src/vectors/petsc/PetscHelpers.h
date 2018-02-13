// This file contains so definitions and wrapper functions for PETSc
#ifndef PETSC_HELPERS
#define PETSC_HELPERS


// Forward declare a few types for PETSc
typedef int PetscErrorCode;
typedef struct _p_Vec* Vec;
typedef struct _p_Mat* Mat;
typedef struct _p_PetscRandom* PetscRandom;



namespace PETSC {


/********************************************************
 * Helper functions                                      *
 ********************************************************/
PetscErrorCode vecDestroy( Vec *v );
PetscErrorCode randomDestroy( PetscRandom *random );
PetscErrorCode matDestroy( Mat *mat );


} // namespace PETSC

#endif
