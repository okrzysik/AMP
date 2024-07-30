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
}

size_t MatrixParameters::getLocalNumberOfRows() const
{
    return d_DOFManagerLeft->numLocalDOF();
}

size_t MatrixParameters::getLocalNumberOfColumns() const
{
    return d_DOFManagerRight->numLocalDOF();
}

size_t MatrixParameters::getGlobalNumberOfRows() const
{
    return d_DOFManagerLeft->numGlobalDOF();
}

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
