
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/boundary/ColumnBoundaryOperator.h"
#include "AMP/operators/boundary/LinearBoundaryOperatorParameters.h"
#include "AMP/utils/Utilities.h"

#include <stdexcept>

namespace AMP {
namespace Operator {

LinearBVPOperator::LinearBVPOperator( const std::shared_ptr<BVPOperatorParameters> &params )
    : LinearOperator( params ),
      d_volumeOperator( std::dynamic_pointer_cast<LinearOperator>( params->d_volumeOperator ) ),
      d_boundaryOperator( params->d_boundaryOperator )
{
    d_Mesh   = d_volumeOperator->getMesh();
    d_matrix = d_volumeOperator->getMatrix();
}

void LinearBVPOperator::reset( const std::shared_ptr<OperatorParameters> &params )
{
    std::shared_ptr<BVPOperatorParameters> inParams =
        std::dynamic_pointer_cast<BVPOperatorParameters>( params );

    AMP_INSIST( ( inParams.get() != nullptr ), "LinearBVPOperator :: reset Null parameter" );

    d_volumeOperator->reset( inParams->d_volumeOperatorParams );

    // first case - single linear boundary operator parameter object
    // This logic does not work with NeumannVectorCorrection boundary
    // operator. As Neumann does not do a matrix correction and its params is
    // not derived from LinearBoundaryOperatorParameters - Allu
    std::shared_ptr<LinearBoundaryOperatorParameters> linearBoundaryParams =
        std::dynamic_pointer_cast<LinearBoundaryOperatorParameters>(
            inParams->d_boundaryOperatorParams );

    if ( linearBoundaryParams != nullptr ) {
        linearBoundaryParams->d_inputMatrix = d_volumeOperator->getMatrix();
        d_boundaryOperator->reset( linearBoundaryParams );
    } else {
        std::shared_ptr<ColumnBoundaryOperatorParameters> columnBoundaryParams =
            std::dynamic_pointer_cast<ColumnBoundaryOperatorParameters>(
                inParams->d_boundaryOperatorParams );

        AMP_ASSERT( columnBoundaryParams != nullptr );

        for ( auto cparams : columnBoundaryParams->d_OperatorParameters ) {

            std::shared_ptr<LinearBoundaryOperatorParameters> linearColBoundaryParams =
                std::dynamic_pointer_cast<LinearBoundaryOperatorParameters>( cparams );
            if ( linearColBoundaryParams != nullptr ) {
                linearColBoundaryParams->d_inputMatrix = d_volumeOperator->getMatrix();
            }
        }
        d_boundaryOperator->reset( columnBoundaryParams );
    }

    d_matrix = d_volumeOperator->getMatrix();
}

void LinearBVPOperator::modifyRHSvector( AMP::LinearAlgebra::Vector::shared_ptr rhs )
{
    ( this->getBoundaryOperator() )->addRHScorrection( rhs );
    ( this->getBoundaryOperator() )->setRHScorrection( rhs );
}
} // namespace Operator
} // namespace AMP
