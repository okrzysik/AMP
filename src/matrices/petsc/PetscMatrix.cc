#include "AMP/matrices/petsc/PetscMatrix.h"
#include "AMP/matrices/petsc/ManagedPetscMatrix.h"


namespace AMP {
namespace LinearAlgebra {


Matrix::shared_ptr PetscMatrix::createView( shared_ptr in_matrix )
{
    auto mat = dynamic_pointer_cast<ManagedPetscMatrix>( in_matrix );
    AMP_INSIST( mat != nullptr, "Managed memory matrix is not well defined" );
    return mat;
}
} // namespace LinearAlgebra
} // namespace AMP
