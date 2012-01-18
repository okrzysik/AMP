
#ifndef included_AMP_ColumnOperator
#define included_AMP_ColumnOperator

#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "vectors/Vector.h"

#include <vector>

namespace AMP {
namespace Operator {

  /**
    A class for representing a composite operator, F=(F1, F2, F3, .., Fk),
    where each of F1,.., Fk are operators themselves. The user is expected to have
    created and initialized the operators F1,.., Fk
    */
  class ColumnOperator : public Operator 
  {

  public :
    // the parameter object for the column operator is intentionally meant not to do
    // anything ColumnOperator specific. Please keep that way
      ColumnOperator(const boost::shared_ptr<OperatorParameters>& params)
        : Operator () { (void) params; }

      /** Default empty constructor */
      ColumnOperator() : Operator() { }

      virtual ~ColumnOperator() { }

      /**
       * The apply routine for the column operator calls apply on each of the component operators
       */
      virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u, 
          AMP::LinearAlgebra::Vector::shared_ptr &r, const double a = -1.0, const double b = 1.0);

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
      virtual void append(boost::shared_ptr< Operator > op);

      /**
       * returns a MultiVariable object corresponding to the ColumnOperator
       * should be called only after all column operators have been appended.
       * no checks to do this right now.
       */
      virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

      virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable();

      bool isValidInput(boost::shared_ptr<AMP::LinearAlgebra::Vector> &u);
      
      boost::shared_ptr< Operator > getOperator(const int i){return d_Operators[i]; }

      int getNumberOfOperators(void){return d_Operators.size(); }
      
    protected :

      std::vector< boost::shared_ptr< Operator > > d_Operators;

    private :

  };

}
}

#endif

