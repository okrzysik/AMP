
#ifndef included_AMP_NeutronicsRhsExtrasParameters
#define included_AMP_NeutronicsRhsExtrasParameters

/* AMP files */
#include "operators/OperatorParameters.h"
#include "utils/Database.h"
#include "vectors/Vector.h"
//#include "ampmesh/MeshUtils.h"

/* Boost files */
#include "boost/shared_ptr.hpp"

namespace AMP {
namespace Operator {

  /**
    A class for encapsulating the parameters that are required for constructing the 
    neutronics source operator. 
    @see NeutronicsRhsExtras
    */
  class NeutronicsRhsExtrasParameters : public OperatorParameters {
    public :
      
      typedef boost::shared_ptr<AMP::Database>  SP_Database;

      NeutronicsRhsExtrasParameters(const SP_Database &db)
    : OperatorParameters(db){} 

//      boost::shared_ptr<AMP::MeshUtils> d_MeshUtils; 

  };

}
}

#endif



