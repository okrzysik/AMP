
#ifndef included_AMP_NeutronicsRhsExtrasParameters
#define included_AMP_NeutronicsRhsExtrasParameters

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
  @see NeutronicsRhsExtras
  */
class NeutronicsRhsExtrasParameters : public OperatorParameters {
public:
    typedef AMP::shared_ptr<AMP::Database> SP_Database;

    explicit NeutronicsRhsExtrasParameters( const SP_Database &db ) : OperatorParameters( db )
    {
        d_numExtras = 0;
    }

    //      AMP::shared_ptr<AMP::MeshUtils> d_MeshUtils;
    int d_numExtras;
    std::vector<std::string> d_extrasName;
};
}
}

#endif
