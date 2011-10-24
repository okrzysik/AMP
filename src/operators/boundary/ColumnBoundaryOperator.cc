
#include "ColumnBoundaryOperator.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Operator {

  void ColumnBoundaryOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr & f,
      const AMP::LinearAlgebra::Vector::shared_ptr & u, AMP::LinearAlgebra::Vector::shared_ptr & r,
      const double a, const double b) {
    for(unsigned int i = 0; i < d_Operators.size(); i++) {
      d_Operators[i]->apply(f, u, r, a, b);
    }
  }

  boost::shared_ptr<OperatorParameters> ColumnBoundaryOperator :: 
    getJacobianParameters(const AMP::LinearAlgebra::Vector::shared_ptr & u) {
      boost::shared_ptr<AMP::Database> db;
      boost::shared_ptr<ColumnBoundaryOperatorParameters> opParameters(new ColumnBoundaryOperatorParameters(db));

      (opParameters->d_OperatorParameters).resize(d_Operators.size());

      for(unsigned int i = 0; i < d_Operators.size(); i++) {
        (opParameters->d_OperatorParameters)[i] = (d_Operators[i]->getJacobianParameters(u));
      }

      return opParameters;
    }

  void ColumnBoundaryOperator :: reset(const boost::shared_ptr<OperatorParameters>& params) {
    boost::shared_ptr<ColumnBoundaryOperatorParameters> columnParameters = 
      boost::dynamic_pointer_cast<ColumnBoundaryOperatorParameters>(params);

    AMP_INSIST( (columnParameters.get() != NULL), "ColumnBoundaryOperator::reset parameter object is NULL" );

    AMP_INSIST( ( ((columnParameters->d_OperatorParameters).size()) == (d_Operators.size()) ),
        " std::vector sizes do not match! " ); 

    for(unsigned int i = 0; i < d_Operators.size(); i++) {
      d_Operators[i]->reset((columnParameters->d_OperatorParameters)[i]);
    }
  }

  void ColumnBoundaryOperator :: append(boost::shared_ptr< BoundaryOperator > op) {
    AMP_INSIST( (op.get() != NULL), "AMP::Operator::ColumnBoundaryOperator::appendRow input argument is a NULL operator");

    d_Operators.push_back(op);
  }

  void ColumnBoundaryOperator :: addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs) {
    for(unsigned int i = 0; i < d_Operators.size(); i++) {
      d_Operators[i]->addRHScorrection(rhs);
    }//end for i
  }

  void ColumnBoundaryOperator :: setRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs) {
    for(unsigned int i = 0; i < d_Operators.size(); i++) {
      d_Operators[i]->setRHScorrection(rhs);
    }//end for i
  }

  void ColumnBoundaryOperator :: modifyInitialSolutionVector(AMP::LinearAlgebra::Vector::shared_ptr sol) {
    for(unsigned int i = 0; i < d_Operators.size(); i++) {
      d_Operators[i]->modifyInitialSolutionVector(sol);
    }//end for i
  }

}
}


