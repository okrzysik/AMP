#include "AMP/matrices/trilinos/EpetraMatrixHelpers.h"
#include "AMP/matrices/Matrix.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/matrices/trilinos/EpetraMatrixData.h"
#include "AMP/matrices/trilinos/ManagedEpetraMatrix.h"


namespace AMP::LinearAlgebra {

/********************************************************
 * Get an Epetra matrix from an AMP matrix              *
 ********************************************************/
std::shared_ptr<ManagedEpetraMatrix> getEpetraMatrix( std::shared_ptr<Matrix> mat )
{
    AMP_ASSERT( mat );
    if ( mat->type() == "ManagedEpetraMatrix" ) {
        return std::dynamic_pointer_cast<ManagedEpetraMatrix>( mat );
    } else {
        // the next piece of code is VERY inefficient and should be optimized in future
        // if we don't deprecate Epetra
        auto matParams = std::make_shared<MatrixParameters>(
            mat->getLeftDOFManager(), mat->getRightDOFManager(), mat->getComm() );

        for ( size_t row = mat->beginRow(); row != mat->endRow(); ++row ) {
            const int i = row - mat->beginRow();
            auto cols   = mat->getColumnIDs( row );
            matParams->setEntriesInRow( i, static_cast<int>( cols.size() ) );
            // matParams->addColumns( cols );
        }

        auto epetraMat = std::make_shared<ManagedEpetraMatrix>( matParams );

        auto data = std::dynamic_pointer_cast<EpetraMatrixData>( epetraMat->getMatrixData() );
        AMP_ASSERT( data );

        std::vector<size_t> cols;
        std::vector<double> vals;
        for ( size_t row = mat->beginRow(); row != mat->endRow(); ++row ) {
            mat->getRowByGlobalID( row, cols, vals );
            data->createValuesByGlobalID( row, cols );
            epetraMat->setValuesByGlobalID( 1, cols.size(), &row, cols.data(), vals.data() );
        }
        epetraMat->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
        return epetraMat;
    }
}


} // namespace AMP::LinearAlgebra
