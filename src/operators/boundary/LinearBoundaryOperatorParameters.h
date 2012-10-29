
#ifndef included_LinearBoundaryOperatorParameters
#define included_LinearBoundaryOperatorParameters

#include "operators/OperatorParameters.h"
#include "matrices/Matrix.h"

namespace AMP {
  namespace Operator {

    class LinearBoundaryOperatorParameters: public OperatorParameters
    {
      public:

        LinearBoundaryOperatorParameters(const boost::shared_ptr<AMP::Database> & db) : 
          OperatorParameters(db){}

        virtual ~LinearBoundaryOperatorParameters(){}

        AMP::LinearAlgebra::Matrix::shared_ptr d_inputMatrix;

    };

  }
}

#endif
