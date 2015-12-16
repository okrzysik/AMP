
#ifndef included_AMP_NeutronicsRhsParameters
#define included_AMP_NeutronicsRhsParameters

/* AMP files */
#include "operators/OperatorParameters.h"
#include "utils/Database.h"
#include "vectors/Vector.h"
//#include "ampmesh/MeshUtils.h"

/* Boost files */
#include "utils/shared_ptr.h"

namespace AMP {
namespace Operator {

/**
  A class for encapsulating the parameters that are required for constructing the
  neutronics source operator.
  @see NeutronicsRhs
  */
class NeutronicsRhsParameters : public OperatorParameters {
public:
    typedef AMP::shared_ptr<AMP::Database> SP_Database;

    explicit NeutronicsRhsParameters( const SP_Database &db ) : OperatorParameters( db ) {}

    //      AMP::shared_ptr<AMP::MeshUtils> d_MeshUtils;
};
}
}

#endif
