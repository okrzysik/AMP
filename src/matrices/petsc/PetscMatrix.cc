#include "matrices/petsc/PetscMatrix.h"
#include "matrices/petsc/ManagedPetscMatrix.h"


namespace AMP {
namespace LinearAlgebra {


Matrix::shared_ptr PetscMatrix::createView( shared_ptr in_matrix )
{
    if ( in_matrix->isA<ManagedPetscMatrix>() ) return in_matrix;

    AMP_ERROR( "Managed memory matrix is not well defined" );
    return Matrix::shared_ptr();
}
}
} // end namespace
