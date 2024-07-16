#include "AMP/matrices/MatrixParameters.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"


namespace AMP::LinearAlgebra {


MatrixParameters::MatrixParameters( std::shared_ptr<AMP::Discretization::DOFManager> left,
                                    std::shared_ptr<AMP::Discretization::DOFManager> right,
                                    const AMP_MPI &comm,
				    const std::function<std::vector<size_t>( size_t )> getRow )
  : MatrixParametersBase( comm ),d_getRowFunction( getRow )
{
    AMP_ASSERT( left );
    AMP_ASSERT( right );
    d_DOFManagerLeft  = left;
    d_DOFManagerRight = right;

    // Create default getRow function if not provided
    // This was migrated to here from createMatrix in MatrixBuilder.cpp
    if ( !d_getRowFunction && d_DOFManagerLeft && d_DOFManagerRight ) {
        d_getRowFunction = [ldof=d_DOFManagerLeft,rdof=d_DOFManagerRight]( size_t row ) {
            auto elem = ldof->getElement( row );
            return rdof->getRowDOFs( elem );
        };
    }
}


size_t MatrixParameters::getLocalNumberOfRows() const { return d_DOFManagerLeft->numLocalDOF(); }


size_t MatrixParameters::getLocalNumberOfColumns() const
{
    return d_DOFManagerRight->numLocalDOF();
}

size_t MatrixParameters::getGlobalNumberOfRows() const { return d_DOFManagerLeft->numGlobalDOF(); }

size_t MatrixParameters::getGlobalNumberOfColumns() const
{
    return d_DOFManagerRight->numGlobalDOF();
}

std::shared_ptr<AMP::Discretization::DOFManager> MatrixParameters::getLeftDOFManager()
{
    return d_DOFManagerLeft;
}

std::shared_ptr<AMP::Discretization::DOFManager> MatrixParameters::getRightDOFManager()
{
    return d_DOFManagerRight;
}

} // namespace AMP::LinearAlgebra
