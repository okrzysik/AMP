#ifndef included_TimeOperatorParameters
#define included_TimeOperatorParameters

#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "utils/shared_ptr.h"

namespace AMP {
namespace TimeIntegrator {

class TimeOperatorParameters : public AMP::Operator::OperatorParameters {
public:
    /**
     * Construct and initialize a parameter list according to input
     * data.  Guess what the required and optional keywords are.
     */
    explicit TimeOperatorParameters( const AMP::shared_ptr<AMP::Database> &db );
    /**
     * Destructor.
     */
    virtual ~TimeOperatorParameters();

    /**
     * Right hand side operator when time operator is written as: u_t = f(u)+g
     * This pointer should be NULL
     * (1) if the parameter object is being used for a reset and not for construction
     */

    AMP::shared_ptr<AMP::Operator::Operator> d_pRhsOperator;

    /**
     * Mass operator which may or may not be present (should be present for FEM)
     * This pointer should be NULL
     * (1) if the parameter object is being used for a reset and not for construction
     */
    AMP::shared_ptr<AMP::Operator::Operator> d_pMassOperator;


    /**
     * Approximate solution at previous time level
     */
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pPreviousTimeSolution;

    /**
     * Source/sink term as well as term containing boundary corrections from mass and rhs operators
     */
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pSourceTerm;

    /**
     * Parameters to reset the rhs operator, this pointer should be NULL only in two cases
     * (1) if we have a linear rhs operator,
     * (2) during construction phase when a non NULL d_pRhsOperator should be supplied
     */
    AMP::shared_ptr<AMP::Operator::OperatorParameters> d_pRhsOperatorParameters;

    /**
     * Parameters to reset the lhs mass operator, this pointer should be NULL only in two cases
     * (1) if we have a linear mass operator,
     * (2) during construction phase when a non NULL d_pMassOperator should be supplied
     */
    AMP::shared_ptr<AMP::Operator::OperatorParameters> d_pMassOperatorParameters;

    /**
     * algebraic variable
     */
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_pAlgebraicVariable;
};
}
}

#endif
