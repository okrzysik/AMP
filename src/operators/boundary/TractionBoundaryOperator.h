
#ifndef included_AMP_TractionBoundaryOperator
#define included_AMP_TractionBoundaryOperator

#include "BoundaryOperator.h"
#include "TractionBoundaryOperatorParameters.h"

namespace AMP {
  namespace Operator {

    class TractionBoundaryOperator : public BoundaryOperator {
      public:

        TractionBoundaryOperator(const boost::shared_ptr<TractionBoundaryOperatorParameters> & params);

        virtual ~TractionBoundaryOperator() { }

        void setTractionVec(AMP::LinearAlgebra::Vector::shared_ptr vec) {
          d_tractionVec = mySubsetVector(vec, d_inputVar);
        }

        void apply(const AMP::LinearAlgebra::Vector::shared_ptr &, const AMP::LinearAlgebra::Vector::shared_ptr &,
            AMP::LinearAlgebra::Vector::shared_ptr &r, const double, const double);

        void addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs);

      protected :

        void computeCorrection();

        AMP::LinearAlgebra::Vector::shared_ptr mySubsetVector(AMP::LinearAlgebra::Vector::shared_ptr vec, 
            AMP::LinearAlgebra::Variable::shared_ptr var);

        AMP::LinearAlgebra::Variable::shared_ptr d_inputVar;
        AMP::LinearAlgebra::Variable::shared_ptr d_outputVar;
        AMP::LinearAlgebra::Vector::shared_ptr d_tractionVec;
        AMP::LinearAlgebra::Vector::shared_ptr d_correction;
        bool d_residualMode;
        short int d_boundaryId;
    };

  }
}

#endif



