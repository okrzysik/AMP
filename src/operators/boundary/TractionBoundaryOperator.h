
#ifndef included_AMP_TractionBoundaryOperator
#define included_AMP_TractionBoundaryOperator

#include "BoundaryOperator.h"
#include "TractionBoundaryOperatorParameters.h"

namespace AMP {
  namespace Operator {

    class TractionBoundaryOperator : public BoundaryOperator {
      public:

        TractionBoundaryOperator(const boost::shared_ptr<TractionBoundaryOperatorParameters> & params) {
          d_variable = params->d_variable;
          d_tractionVec = mySubsetVector(params->d_tractionVec, params->d_variable);
          reset(params);
        }

        virtual ~TractionBoundaryOperator() { }

        void setVariable(const AMP::LinearAlgebra::Variable::shared_ptr & var) {
          d_variable = var;
        }

        void setTractionVec(AMP::LinearAlgebra::Vector::shared_ptr vec) {
          d_tractionVec = mySubsetVector(vec, d_variable);
        }

        virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &, const AMP::LinearAlgebra::Vector::shared_ptr &,
            AMP::LinearAlgebra::Vector::shared_ptr &r, const double, const double) {
          if(d_residualMode) {
            AMP::LinearAlgebra::Vector::shared_ptr rInternal = mySubsetVector(r, d_variable);
            rInternal->subtract(rInternal, d_correction);
          }
        }

        virtual void addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs) {
          if(!d_residualMode) {
            AMP::LinearAlgebra::Vector::shared_ptr myRhs = mySubsetVector(rhs, d_variable);
            myRhs->add(myRhs, d_correction);
          }
        }

      protected :

        AMP::LinearAlgebra::Vector::shared_ptr mySubsetVector(AMP::LinearAlgebra::Vector::shared_ptr vec, 
            AMP::LinearAlgebra::Variable::shared_ptr var) {
          if(d_Mesh.get() != NULL) {
            AMP::LinearAlgebra::VS_Mesh meshSelector(var->getName(), d_Mesh);
            AMP::LinearAlgebra::Vector::shared_ptr meshSubsetVec = vec->select(meshSelector, var->getName());
            return meshSubsetVec->subsetVectorForVariable(var);
          } else {
            return vec->subsetVectorForVariable(var);
          }
        }

        AMP::LinearAlgebra::Vector::shared_ptr d_tractionVec;

        AMP::LinearAlgebra::Vector::shared_ptr d_correction;

        bool d_residualMode;

        short int d_boundaryId;
    };

  }
}

#endif



