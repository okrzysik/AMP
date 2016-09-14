namespace AMP {
namespace LinearAlgebra {



inline PetscVector::PetscVector() { d_PetscRandom = 0; }


inline PetscRandom &PetscVector::getPetscRandom( AMP_MPI comm )
{
    if ( d_PetscRandom == 0 ) {
        d_PetscRandom = new PetscRandom;
        PetscRandomCreate( comm.getCommunicator(), d_PetscRandom );
        PetscRandomSetType( *d_PetscRandom, PETSCRAND48 ); // This is a horrible RNG for
                                                           // stochastic simulation.  Do not
                                                           // use.
    }
    return *d_PetscRandom;
}


inline Vec &PetscVector::getVec() { return d_petscVec; }


inline Vec PetscVector::getVec() const { return d_petscVec; }


inline PetscVector::~PetscVector()
{
    if ( d_PetscRandom ) {
        PetscVector::RandomDestroy( d_PetscRandom );
        delete d_PetscRandom;
    }
}


/********************************************************
* Helper functions                                      *
********************************************************/
inline PetscErrorCode PetscVector::VecDestroy( Vec *v )
{
    #if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
        return ::VecDestroy( *v );
    #elif PETSC_VERSION_GE(3,2,0)
        return ::VecDestroy( v );
    #else
        #error Not programmed for this version yet
    #endif
}
inline PetscErrorCode PetscVector::RandomDestroy( PetscRandom *random )
{
    #if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
        return ::PetscRandomDestroy( *random );
    #elif PETSC_VERSION_GE(3,2,0)
        return ::PetscRandomDestroy( random );
    #else
        #error Not programmed for this version of petsc
    #endif
}


} // LinearAlgebra
} // AMP

