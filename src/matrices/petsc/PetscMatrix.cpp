#include "AMP/matrices/petsc/PetscMatrix.h"
#include "AMP/matrices/petsc/ManagedPetscMatrix.h"


namespace AMP {
namespace LinearAlgebra {


Matrix::shared_ptr PetscMatrix::createView( shared_ptr in_matrix )
{
#ifdef USE_EXT_TRILINOS
    auto mat = std::dynamic_pointer_cast<ManagedPetscMatrix>( in_matrix );
    AMP_INSIST( mat != nullptr, "Managed memory matrix is not well defined" );
    return mat;
#else
    NULL_USE( in_matrix );
    AMP_ERROR( "ManagedPetscMatrix currently requires Epetra" );
    return nullptr;
#endif
}
} // namespace LinearAlgebra
} // namespace AMP
