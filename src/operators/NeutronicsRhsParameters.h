#ifndef included_AMP_NeutronicsRhsParameters
#define included_AMP_NeutronicsRhsParameters

#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/Database.h"
#include "AMP/vectors/Vector.h"
#include <memory>

namespace AMP {
namespace Operator {

/**
  A class for encapsulating the parameters that are required for constructing the
  neutronics source operator.
  @see NeutronicsRhs
  */
class NeutronicsRhsParameters : public OperatorParameters
{
public:
    explicit NeutronicsRhsParameters( std::shared_ptr<AMP::Database> db ) : OperatorParameters( db )
    {
    }
};

} // namespace Operator
} // namespace AMP

#endif
