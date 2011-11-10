
#ifndef included_AMP_LinearOperator
#define included_AMP_LinearOperator

#include "boost/shared_ptr.hpp"
#include "matrices/Matrix.h"
#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "vectors/Vector.h"

namespace AMP {
  namespace Operator {

    /**
      An abstract base class for representing a linear operator. This class 
      stores the matrix representation of the linear operator. It provides
      an implementation of the apply() function.
      @see Operator
      */
    class LinearOperator : public Operator 
    {

      public :

        /**
          Constructor. This resets the matrix shared pointer.
          @param [in] params 
          */
        LinearOperator (const boost::shared_ptr<OperatorParameters> & params);

        /**
          Destructor
          */
        virtual ~LinearOperator() { }

        /**
          The apply function for this operator, A, performs the following operation:
          r = a*A(u) + b*f, if f is not NULL and r = a*A(u), if f is NULL.
          Here, A(u) is simply a Matrix-Vector multiplication.
          @param [in] f auxillary/rhs vector. 
          @param [in] u input vector. 
          @param [out] r residual/output vector. 
          @param [in] a first constant used in the expression: r = a*A(u) + b*f. The default value is -1.
          @param [in] b second constant used in the expression: r = a*A(u) + b*f. The default value is 1.
          */
        virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

        /**
          @return The matrix representation of this linear operator.
          */
        boost::shared_ptr<AMP::LinearAlgebra::Matrix> getMatrix();

        /**
          Copies the shared pointer for the matrix representation of this linear operator.
          @param [in] in_mat The matrix representation of this linear operator.
          */
        virtual void setMatrix(const boost::shared_ptr<AMP::LinearAlgebra::Matrix> & in_mat);

        virtual void resetApplyCount()  {
          d_applyCount = 0;      
        }

        virtual unsigned int getApplyCount() {
          return d_applyCount;
        }

      protected :

        boost::shared_ptr<AMP::LinearAlgebra::Matrix> d_matrix; /**< The matrix shared pointer. */

        unsigned int d_applyCount;

      private :

    };

  }
}

#endif


