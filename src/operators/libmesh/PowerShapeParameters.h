
#ifndef included_AMP_PowerShapeParameters
#define included_AMP_PowerShapeParameters

/* AMP files */
#include "operators/OperatorParameters.h"
#include "utils/Database.h"
#include "vectors/Vector.h"

/* Boost files */
#include "boost/shared_ptr.hpp"

namespace AMP {
namespace Operator {

  /**
    A class for encapsulating the parameters that are required for constructing the 
    PowerShape operator. 
    @see PowerShape
    */
  class PowerShapeParameters : public OperatorParameters {
    public :
      
      PowerShapeParameters(const boost::shared_ptr<AMP::Database> &db)
    : OperatorParameters(db){} 

  };

}
}

#endif



