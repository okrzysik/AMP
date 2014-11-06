
#ifndef included_AMP_OperatorWrapper
#define included_AMP_OperatorWrapper

#include "operators/Operator.h"
#include "vectors/MultiVector.h"
#include "utils/Utilities.h"

namespace AMP {
  namespace Operator {

    class OperatorWrapper : public Operator {
      public:
        OperatorWrapper() { } 

        virtual ~OperatorWrapper() { }

        void setOperator(Operator::shared_ptr op) {
          d_operator = op;
        }

        void reset(const AMP::shared_ptr<OperatorParameters>& params) {
          d_operator->reset(params);
        }

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
          return (d_operator->getOutputVariable());
        }

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1) {
          return (d_operator->getInputVariable(varId));
        }

        AMP::Mesh::MeshManager::Adapter::shared_ptr getMeshAdapter() {
          return (d_operator->getMeshAdapter());
        }

        AMP::shared_ptr<OperatorParameters> 
          getJacobianParameters(const AMP::shared_ptr<AMP::LinearAlgebra::Vector>& u) {
            AMP::LinearAlgebra::Vector::shared_ptr uTmp = createCorrectInputVector(u);
            return (d_operator->getJacobianParameters(uTmp));
          }

        bool isValidInput(AMP::shared_ptr<AMP::LinearAlgebra::Vector>& u) {
          AMP::LinearAlgebra::Vector::shared_ptr uTmp = createCorrectInputVector(u);
          return (d_operator->isValidInput(uTmp));
        }

        void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, 
            AMP::LinearAlgebra::Vector::const_shared_ptr u, AMP::LinearAlgebra::Vector::shared_ptr r,
            const double a = -1.0, const double b = 1.0) {
          AMP::LinearAlgebra::Vector::shared_ptr uTmp = createCorrectInputVector(u);
          d_operator->apply(f, uTmp, r, a, b);
        }

        void setFullVector(AMP::Vector::shared_ptr u) {
          d_fullVector = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>(
              AMP::LinearAlgebra::MultiVector::view(u, u->getComm()));
        }

      protected:
        AMP::LinearAlgebra::Vector::shared_ptr createCorrectInputVector(AMP::LinearAlgebra::Vector::shared_ptr u) {
          AMP_ASSERT(d_fullVector != NULL);
          AMP::Vector::shared_ptr oldVec = d_fullVector->subsetVectorForVariable(getOutputVariable());
          AMP::Vector::shared_ptr newVec = u->subsetVectorForVariable(getOutputVariable());
          d_fullVector->replaceSubVector(oldVec, newVec);
          return d_fullVector;
        }

        Operator::shared_ptr d_operator;
        AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> d_fullVector;
    };

  }
}

#endif


