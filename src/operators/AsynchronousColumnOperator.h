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
      AsynchronousColumnOperator ( const boost::shared_ptr < OperatorParameters > & );

      virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u, 
          AMP::LinearAlgebra::Vector::shared_ptr &r, const double a = -1.0, const double b = 1.0);

      virtual void applyFinish(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u, 
             AMP::LinearAlgebra::Vector::shared_ptr &r, const double a = -1.0, const double b = 1.0);

      virtual void applyStart(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u, 
             AMP::LinearAlgebra::Vector::shared_ptr &r, const double a = -1.0, const double b = 1.0);

      virtual void append(boost::shared_ptr< Operator > op);
  };

}
}

#endif
