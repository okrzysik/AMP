#ifndef included_AMP_BoundaryOperator
#define included_AMP_BoundaryOperator

#include "AMP/operators/Operator.h"

namespace AMP {
namespace Operator {


//  An abstract base class for representing a linear operator.
class BoundaryOperator : public Operator
{

public:
    explicit BoundaryOperator( const std::shared_ptr<OperatorParameters> &params )
        : Operator( params )
    {
    }

    virtual ~BoundaryOperator() {}

    virtual void addRHScorrection( AMP::LinearAlgebra::Vector::shared_ptr ) {}

    virtual void setRHScorrection( AMP::LinearAlgebra::Vector::shared_ptr ) {}

    virtual void modifyInitialSolutionVector( AMP::LinearAlgebra::Vector::shared_ptr ) {}

protected:
private:
};


} // Namespace Operator
} // Namespace AMP

#endif
