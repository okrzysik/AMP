
namespace AMP {
namespace LinearAlgebra {

inline PetscMatrix::PetscMatrix() {}

inline PetscMatrix::PetscMatrix( MatrixParameters::shared_ptr params ) : Matrix( params ) {}

inline PetscMatrix::~PetscMatrix() {}

inline Mat &PetscMatrix::getMat() { return d_Mat; }

inline Mat PetscMatrix::getMat() const { return d_Mat; }


/********************************************************
* Helper functions                                      *
********************************************************/
inline PetscErrorCode PetscMatrix::MatDestroy( Mat *mat )
{
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    return MatDestroy( mat );
#elif PETSC_VERSION_GE(3,2,0)
    return MatDestroy( mat );
#else
    #error Not programmed for this version yet
#endif
}

} // LinearAlgebra
} // AMP
