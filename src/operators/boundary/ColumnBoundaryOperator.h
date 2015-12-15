
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
      explicit ColumnBoundaryOperator(const AMP::shared_ptr<OperatorParameters>& params)
        : BoundaryOperator (params) { }

      virtual ~ColumnBoundaryOperator() { }

      /**
       * The apply routine for the column operator calls apply on each of the component operators
       */
      virtual void apply(AMP::LinearAlgebra::Vector::const_shared_ptr u,
			 AMP::LinearAlgebra::Vector::shared_ptr r) override;

      /**
        A function for computing the information necessary to construct the jacobian.
        @param u The solution vector that is used to construct the jacobian
        @return The parameters required to construct the jacobian.
        */
      AMP::shared_ptr<OperatorParameters> getParameters(const std::string &type,
                                                        AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                                        AMP::shared_ptr<OperatorParameters> params = NULL) override;

      virtual void reset(const AMP::shared_ptr<OperatorParameters>& params);

      /**
       * \param op
       *            shared pointer to an operator to append to the existing column of operators
       */
      virtual void append(AMP::shared_ptr< BoundaryOperator > op);

      AMP::shared_ptr< BoundaryOperator > getBoundaryOperator(int i){ return d_Operators[i]; }

      size_t  numberOfBoundaryOperators () { return d_Operators.size(); }

      void addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr );

      void setRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr );

      void modifyInitialSolutionVector(AMP::LinearAlgebra::Vector::shared_ptr );

    protected :

      std::vector< AMP::shared_ptr< BoundaryOperator > > d_Operators;

    private :

  };

}
}

#endif

