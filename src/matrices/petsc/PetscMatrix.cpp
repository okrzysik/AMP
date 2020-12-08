#include "AMP/matrices/petsc/PetscMatrix.h"

#include "petsc/private/petscimpl.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * view                                                          *
 ****************************************************************/
std::shared_ptr<const PetscMatrix> PetscMatrix::constView( Matrix::const_shared_ptr inMat )
{
    return view( std::const_pointer_cast<Matrix>( inMat ) );
}
std::shared_ptr<PetscMatrix> PetscMatrix::view( Matrix::shared_ptr inMat )
{
    // Check if we have an existing view
    if ( std::dynamic_pointer_cast<PetscMatrix>( inMat ) )
        return std::dynamic_pointer_cast<PetscMatrix>( inMat );
    /*if ( inMat->hasView<PetscMatrix>() )
        return inMat->getView<PetscMatrix>();*/
    // Create the view
    std::shared_ptr<PetscMatrix> ptr( new PetscMatrix( inMat ) );
    // inMat->registerView( ptr );
    return ptr;
}

PetscMatrix::PetscMatrix() {}
PetscMatrix::PetscMatrix( std::shared_ptr<Matrix> mat )
    : d_Mat( PETSC::getMat( mat ) ), d_matrix( mat )
{
}
PetscMatrix::~PetscMatrix() { PETSC::matDestroy( &d_Mat ); }


} // namespace LinearAlgebra
} // namespace AMP
