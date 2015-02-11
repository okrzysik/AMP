
#include "operators/ConstraintsEliminationOperator.h"

#include <algorithm>
#include <cmath>

namespace AMP {
  namespace Operator {

    void ConstraintsEliminationOperator::reset(const AMP::shared_ptr<OperatorParameters> & params) {

      AMP_INSIST( (params != NULL), "NULL parameter" );
      AMP_INSIST( ((params->d_db) != NULL), "NULL database" );

    }

    AMP::LinearAlgebra::Variable::shared_ptr ConstraintsEliminationOperator::getInputVariable() {
      return d_InputVariable;
    }

    AMP::LinearAlgebra::Variable::shared_ptr ConstraintsEliminationOperator::getOutputVariable() {
      return d_OutputVariable;
    }

    void ConstraintsEliminationOperator::apply(AMP::LinearAlgebra::Vector::const_shared_ptr, 
        AMP::LinearAlgebra::Vector::const_shared_ptr,
        AMP::LinearAlgebra::Vector::shared_ptr r, const double, const double ) {
      addSlaveToMaster(r);
      setSlaveToZero(r);
    }

    void ConstraintsEliminationOperator::setSlaveToZero(AMP::LinearAlgebra::Vector::shared_ptr u) {
      if (!d_SlaveIndices.empty()) {
        std::vector<double> zeroSlaveValues(d_SlaveIndices.size(), 0.0);
        u->setLocalValuesByGlobalID(d_SlaveIndices.size(), &(d_SlaveIndices[0]), &(zeroSlaveValues[0]));
      } // end if
    }

    void ConstraintsEliminationOperator::addShiftToSlave(AMP::LinearAlgebra::Vector::shared_ptr u) {
      AMP_ASSERT( d_SlaveShift.size() == d_SlaveIndices.size() );
      if (!d_SlaveIndices.empty()) {
        u->addLocalValuesByGlobalID(d_SlaveIndices.size(), &(d_SlaveIndices[0]), &(d_SlaveShift[0])); 
      } // end if
    }

  } // end namespace Operator
} // end namespace AMP


