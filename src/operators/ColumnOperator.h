
#ifndef included_AMP_ColumnOperator
#define included_AMP_ColumnOperator

#include "AMP/operators/ColumnOperatorParameters.h"
#include "AMP/operators/Operator.h"

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

public:
    // the parameter object for the column operator is intentionally meant not to do
    // anything ColumnOperator specific. Please keep that way
    explicit ColumnOperator( const std::shared_ptr<OperatorParameters> &params ) : Operator()
    {
        (void) params;
    }

    /** Default empty constructor */
    ColumnOperator() : Operator() {}

    virtual ~ColumnOperator() {}

    /**
     * The apply routine for the column operator calls apply on each of the component operators
     */
    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr r ) override;

    /**
     * Column specific implementation of the residual: f-L(u)
     * \param f: shared pointer to const vector rhs
     * \param u: shared pointer to const vector u
     * \param r: shared pointer to vector residual
     */
    virtual void residual( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                           AMP::LinearAlgebra::Vector::const_shared_ptr u,
                           AMP::LinearAlgebra::Vector::shared_ptr r ) override;

    virtual void reset( const std::shared_ptr<OperatorParameters> &params ) override;

    /**
      A function for computing the information necessary to construct the jacobian.
     * \param type: std:string specifying type of return operator parameters
     *      being requested. Currently the valid option is Jacobian
     * \param u: const pointer to current solution vector
     * \param params: pointer to additional parameters that might be required
     *      to construct the return parameters
      @return The parameters required to construct the jacobian.
      */
    virtual std::shared_ptr<OperatorParameters>
    getParameters( const std::string &type,
                   AMP::LinearAlgebra::Vector::const_shared_ptr u,
                   std::shared_ptr<OperatorParameters> params = nullptr ) override;

    /**
     * \param op
     *            shared pointer to an operator to append to the existing column of operators
     */
    virtual void append( std::shared_ptr<Operator> op );

    /**
     * returns a MultiVariable object corresponding to the ColumnOperator
     * should be called only after all column operators have been appended.
     * no checks to do this right now.
     */
    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override;

    virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override;

    bool isValidInput( std::shared_ptr<AMP::LinearAlgebra::Vector> &u ) override;

    std::shared_ptr<Operator> getOperator( size_t i ) { return d_Operators[i]; }

    size_t getNumberOfOperators( void ) { return d_Operators.size(); }

protected:
    std::vector<std::shared_ptr<Operator>> d_Operators;

private:
};
} // namespace Operator
} // namespace AMP

#endif
