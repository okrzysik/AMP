
#ifndef included_AMP_BoundaryOperator
#define included_AMP_BoundaryOperator

#include "operators/Operator.h"

namespace AMP {
namespace Operator {

  //  An abstract base class for representing a linear operator.
  class BoundaryOperator : public Operator 
  {

    public :

      BoundaryOperator (const boost::shared_ptr<OperatorParameters> & params)
        : Operator (params) { }

      virtual ~BoundaryOperator() { }

      virtual void addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr ) {
        //Do nothing
      }

      virtual void setRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr ) {
        //Do nothing
      }

      virtual void modifyInitialSolutionVector(AMP::LinearAlgebra::Vector::shared_ptr ) {
        //Do nothing
      }

    protected :

    private :

  };

}
}

#endif


