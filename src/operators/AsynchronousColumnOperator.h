#ifndef included_AMP_AsynchronousColumnOperator
#define included_AMP_AsynchronousColumnOperator

#include "ColumnOperator.h"

namespace AMP {
namespace Operator {

  /** \brief  A column operator of asynchronous operators.  The apply method will start the list
      of operators then finalize the list of operators
      */
  class AsynchronousColumnOperator : public ColumnOperator
  {
    public:
      /** Constructor
        */
      AsynchronousColumnOperator ( const AMP::shared_ptr < OperatorParameters > & );

      virtual void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u, 
          AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0);

      virtual void applyFinish(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u, 
             AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0);

      virtual void applyStart(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u, 
             AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0);

      virtual void append(AMP::shared_ptr< Operator > op);
  };

}
}

#endif
