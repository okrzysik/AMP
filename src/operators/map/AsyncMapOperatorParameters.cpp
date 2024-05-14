#include "AMP/operators/map/AsyncMapOperatorParameters.h"

namespace AMP::Operator {


AsyncMapOperatorParameters::AsyncMapOperatorParameters( std::shared_ptr<AMP::Database> db )
    : AsynchronousOperatorParameters( db )
{
    d_BoundaryID1         = -1;
    d_BoundaryID2         = -1;
    d_commTag             = -1;
    callMakeConsistentSet = true;
}


AsyncMapOperatorParameters::~AsyncMapOperatorParameters() = default;
} // namespace AMP::Operator
