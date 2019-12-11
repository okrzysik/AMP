#ifndef included_AMP_PowerShapeParameters
#define included_AMP_PowerShapeParameters

#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Vector.h"

namespace AMP {
namespace Operator {

/**
  A class for encapsulating the parameters that are required for constructing the
  PowerShape operator.
  @see PowerShape
  */
class PowerShapeParameters : public OperatorParameters
{
public:
    explicit PowerShapeParameters( AMP::shared_ptr<AMP::Database> db ) : OperatorParameters( db ) {}
};
} // namespace Operator
} // namespace AMP

#endif
