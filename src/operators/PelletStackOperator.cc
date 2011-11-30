
#include "operators/PelletStackOperator.h"

namespace AMP {
  namespace Operator {

    PelletStackOperator :: PelletStackOperator(const boost::shared_ptr<OperatorParameters> & params)
      : Operator(params) {
      }

    void PelletStackOperator :: setCurrentPellet(unsigned int pellId) {
    }

    unsigned int PelletStackOperator :: getTotalNumberOfPellets() {
      return 0;
    }

    bool PelletStackOperator :: hasPellet(unsigned int pellId) {
      return true;
    }

    void PelletStackOperator :: applyScaling(AMP::LinearAlgebra::Vector::shared_ptr f) {
    }

    void PelletStackOperator :: applyUnscaling(AMP::LinearAlgebra::Vector::shared_ptr f) {
    }

    AMP::LinearAlgebra::Variable::shared_ptr PelletStackOperator :: getOutputVariable() {
      AMP::LinearAlgebra::Variable::shared_ptr emptyPointer;
      return emptyPointer;
    }

    AMP::LinearAlgebra::Variable::shared_ptr PelletStackOperator :: getInputVariable(int varId) {
      AMP::LinearAlgebra::Variable::shared_ptr emptyPointer;
      return emptyPointer;
    }

    void PelletStackOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
        const double a, const double b) {
    }

  }
}


