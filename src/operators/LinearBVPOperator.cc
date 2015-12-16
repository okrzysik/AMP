
#include "operators/LinearBVPOperator.h"
#include "operators/boundary/ColumnBoundaryOperator.h"
#include "operators/boundary/LinearBoundaryOperatorParameters.h"
#include "utils/Utilities.h"

#include <stdexcept>

namespace AMP {
namespace Operator {

LinearBVPOperator::LinearBVPOperator( const AMP::shared_ptr<BVPOperatorParameters> &params )
    : LinearOperator( params )
{
    d_volumeOperator   = AMP::dynamic_pointer_cast<LinearOperator>( params->d_volumeOperator );
    d_boundaryOperator = params->d_boundaryOperator;
    d_Mesh             = d_volumeOperator->getMesh();
    d_matrix           = d_volumeOperator->getMatrix();
}

void LinearBVPOperator::reset( const AMP::shared_ptr<OperatorParameters> &params )
{
    AMP::shared_ptr<BVPOperatorParameters> inParams =
        AMP::dynamic_pointer_cast<BVPOperatorParameters>( params );

    AMP_INSIST( ( inParams.get() != NULL ), "LinearBVPOperator :: reset Null parameter" );

    d_volumeOperator->reset( inParams->d_volumeOperatorParams );

    // first case - single linear boundary operator parameter object
    // This logic does not work with NeumannVectorCorrection boundary
    // operator. As Neumann does not do a matrix correction and its params is
    // not derived from LinearBoundaryOperatorParameters - Allu
    AMP::shared_ptr<LinearBoundaryOperatorParameters> linearBoundaryParams =
        AMP::dynamic_pointer_cast<LinearBoundaryOperatorParameters>(
            inParams->d_boundaryOperatorParams );

    if ( linearBoundaryParams != NULL ) {
        linearBoundaryParams->d_inputMatrix = d_volumeOperator->getMatrix();
        d_boundaryOperator->reset( linearBoundaryParams );
    }
    else {
        AMP::shared_ptr<ColumnBoundaryOperatorParameters> columnBoundaryParams =
            AMP::dynamic_pointer_cast<ColumnBoundaryOperatorParameters>(
                inParams->d_boundaryOperatorParams );

        AMP_ASSERT( columnBoundaryParams != NULL );

        for ( unsigned int i = 0; i < columnBoundaryParams->d_OperatorParameters.size(); i++ ) {
            AMP::shared_ptr<OperatorParameters> cparams =
                columnBoundaryParams->d_OperatorParameters[i];
            AMP::shared_ptr<LinearBoundaryOperatorParameters> linearBoundaryParams =
                AMP::dynamic_pointer_cast<LinearBoundaryOperatorParameters>( cparams );
            if ( linearBoundaryParams != NULL ) {
                linearBoundaryParams->d_inputMatrix = d_volumeOperator->getMatrix();
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
}
}
