#ifndef included_AMP_BVPOperatorParameters
#define included_AMP_BVPOperatorParameters

#include "AMP/operators/OperatorParameters.h"
#include "AMP/operators/boundary/BoundaryOperator.h"
#include <memory>


namespace AMP::Operator {


/** \class BVPOperatorParameters
 *  \brief Parameter object used for BVPOperators both linear and nonlinear
 *  \details This parameter object is used in the initialization and reset
 *           of NonlinearBVPOperator and LinearBVPOperator
 */

class BVPOperatorParameters : public OperatorParameters
{
public:
    explicit BVPOperatorParameters( std::shared_ptr<AMP::Database> db ) : OperatorParameters( db )
    {
    }

    virtual ~BVPOperatorParameters() {}

    std::shared_ptr<OperatorParameters> d_volumeOperatorParams;

    std::shared_ptr<OperatorParameters> d_boundaryOperatorParams;

    std::shared_ptr<Operator> d_volumeOperator;

    std::shared_ptr<BoundaryOperator> d_boundaryOperator;
};


} // namespace AMP::Operator

#endif
