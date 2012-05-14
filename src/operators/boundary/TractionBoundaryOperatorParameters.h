
#ifndef included_AMP_TractionBoundaryOperatorParameters
#define included_AMP_TractionBoundaryOperatorParameters

#include "operators/OperatorParameters.h"

namespace AMP {
  namespace Operator {

    class TractionBoundaryOperatorParameters : public OperatorParameters {
      public :

        TractionBoundaryOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
          : OperatorParameters(db) {  }

        ~TractionBoundaryOperatorParameters() { }

        AMP::LinearAlgebra::Vector::shared_ptr d_tractionVec;

    };

  }
}

#endif

