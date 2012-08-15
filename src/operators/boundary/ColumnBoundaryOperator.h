
#ifndef included_AMP_ColumnBoundaryOperator
#define included_AMP_ColumnBoundaryOperator

#include "BoundaryOperator.h"
#include "operators/ColumnOperatorParameters.h"
#include "vectors/Vector.h"

#include <vector>

namespace AMP {
namespace Operator {

  /**
    A class for representing a composite boundary operator, F=(F1, F2, F3, .., Fk),
    where each of F1,.., Fk are boundary operators themselves. The user is expected to have
    created and initialized the operators F1,.., Fk. This class can be used to impose a mix of boundary
    conditions for a volume operator over the boundary
    */

  typedef ColumnOperatorParameters ColumnBoundaryOperatorParameters;

  class ColumnBoundaryOperator : public BoundaryOperator
  {

    public :
      ColumnBoundaryOperator(const boost::shared_ptr<OperatorParameters>& params)
        : BoundaryOperator (params) { }

      virtual ~ColumnBoundaryOperator() { }

      /**
       * The apply routine for the column operator calls apply on each of the component operators
       */
      virtual void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
          AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0);

      /**
        A function for computing the information necessary to construct the jacobian.
        @param u The solution vector that is used to construct the jacobian
        @return The parameters required to construct the jacobian.
        */
      virtual boost::shared_ptr<OperatorParameters> getJacobianParameters(const AMP::LinearAlgebra::Vector::shared_ptr & u);

      virtual void reset(const boost::shared_ptr<OperatorParameters>& params);

      /**
       * \param op
       *            shared pointer to an operator to append to the existing column of operators
       */
      virtual void append(boost::shared_ptr< BoundaryOperator > op);

      boost::shared_ptr< BoundaryOperator > getBoundaryOperator(int i){ return d_Operators[i]; }

      size_t  numberOfBoundaryOperators () { return d_Operators.size(); }

      void addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr );

      void setRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr );

      void modifyInitialSolutionVector(AMP::LinearAlgebra::Vector::shared_ptr );

    protected :

      std::vector< boost::shared_ptr< BoundaryOperator > > d_Operators;

    private :

  };

}
}

#endif

