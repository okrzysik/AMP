#ifndef included_AMP_BVPOperatorParameters
#define included_AMP_BVPOperatorParameters

#include "AMP/operators/OperatorParameters.h"
#include "AMP/operators/boundary/BoundaryOperator.h"
#include "AMP/utils/shared_ptr.h"


namespace AMP {
namespace Operator {


/** \class BVPOperatorParameters
 *  \brief Parameter object used for BVPOperators both linear and nonlinear
 *  \details This parameter object is used in the initialization and reset
 *           of NonlinearBVPOperator and LinearBVPOperator
 */

class BVPOperatorParameters : public OperatorParameters
{
public:
    explicit BVPOperatorParameters( const AMP::shared_ptr<AMP::Database> &db )
        : OperatorParameters( db )
    {
    }

    virtual ~BVPOperatorParameters() {}

    AMP::shared_ptr<OperatorParameters> d_volumeOperatorParams;

    AMP::shared_ptr<OperatorParameters> d_boundaryOperatorParams;

    AMP::shared_ptr<Operator> d_volumeOperator;

    AMP::shared_ptr<BoundaryOperator> d_boundaryOperator;
};


} // namespace Operator
} // namespace AMP

#endif
