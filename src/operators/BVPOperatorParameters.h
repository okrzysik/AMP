
#ifndef included_AMP_BVPOperatorParameters
#define included_AMP_BVPOperatorParameters

/*AMP files */
#include "operators/OperatorParameters.h"
#include "LinearOperator.h"
#include "operators/boundary/BoundaryOperator.h"

/*Boost files */
#include "boost/shared_ptr.hpp"

namespace AMP {
namespace Operator {
  /** \class BVPOperatorParameters
   *  \brief Parameter object used for BVPOperators both linear and nonlinear
   *  \details This parameter object is used in the initialization and reset
   *           of NonlinearBVPOperator and LinearBVPOperator
   */
  
  class BVPOperatorParameters : public OperatorParameters {
  public :

      BVPOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : OperatorParameters(db) { }

      virtual ~BVPOperatorParameters() { }

      boost::shared_ptr<OperatorParameters> d_volumeOperatorParams;

      boost::shared_ptr<OperatorParameters> d_boundaryOperatorParams;

      boost::shared_ptr<Operator> d_volumeOperator;

      boost::shared_ptr<BoundaryOperator> d_boundaryOperator;

  };

}  
}

#endif



