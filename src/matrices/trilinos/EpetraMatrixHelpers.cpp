#include "AMP/matrices/trilinos/EpetraMatrixHelpers.h"
#include "AMP/matrices/Matrix.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/trilinos/ManagedEpetraMatrix.h"


namespace AMP::LinearAlgebra {

/********************************************************
 * Get an Epetra vector from an AMP vector               *
 ********************************************************/
std::shared_ptr<ManagedEpetraMatrix> getEpetraMatrix( std::shared_ptr<Matrix> mat )
{
    AMP_ASSERT( mat );
    if ( mat->type() == "ManagedEpetraMatrix" ) {
        return std::dynamic_pointer_cast<ManagedEpetraMatrix>( mat );
    } else {
        auto matParams = std::make_shared<MatrixParameters>(
            mat->getLeftDOFManager(), mat->getRightDOFManager(), mat->getComm() );
        AMP_ERROR( "Not implemented" );
        return std::make_shared<ManagedEpetraMatrix>( matParams );
    }
}


} // namespace AMP::LinearAlgebra
