#ifndef included_AMP_LinearTimeOperator
#define included_AMP_LinearTimeOperator

#include "AMP/operators/LinearOperator.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Vector.h"


namespace AMP {
namespace TimeIntegrator {

/*!
  @brief base class for operator class associated with ImplicitTimeIntegrator

  Class ImplicitLinearTimeOperator is a base class derived from Operator. It
  is the operator class associated with a ImplicitTimeIntegrator. The solver associated
  with the ImplicitTimeIntegrator will register this object.

  @see ImplicitTimeIntegrator
  @see Operator
  @see SolverStrategy
*/
class LinearTimeOperator : public AMP::Operator::LinearOperator
{
public:
    explicit LinearTimeOperator( AMP::shared_ptr<AMP::Operator::OperatorParameters> params );
    virtual ~LinearTimeOperator();

    /**
     * This function is useful for re-initializing an operator
     * \param params
     *        parameter object containing parameters to change
     */
    virtual void reset( const AMP::shared_ptr<AMP::Operator::OperatorParameters> &params ) override;

    void registerRhsOperator( AMP::shared_ptr<AMP::Operator::LinearOperator> op )
    {
        d_pRhsOperator = op;
    }
    void registerMassOperator( AMP::shared_ptr<AMP::Operator::LinearOperator> op )
    {
        d_pMassOperator = op;
    }

    AMP::shared_ptr<Operator> getRhsOperator( void ) { return d_pRhsOperator; }
    AMP::shared_ptr<Operator> getMassOperator( void ) { return d_pMassOperator; }

    void setPreviousSolution( AMP::shared_ptr<AMP::LinearAlgebra::Vector> previousSolution )
    {
        d_pPreviousTimeSolution = previousSolution;
    }

    void setDt( double dt ) { d_dCurrentDt = dt; }

    // added by JL
    void setScalingFactor( double scalingFactor ) { d_dScalingFactor = scalingFactor; }

    // added by JL //correction by RS
    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override
    {
        return d_pRhsOperator->getInputVariable();
    }

    /**
     * returns a Variable object corresponding to the rhs operator
     */
    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override
    {
        return d_pRhsOperator->getOutputVariable();
    }
    // JL
    void registerCurrentTime( double currentTime ) { d_current_time = currentTime; }

    AMP::shared_ptr<AMP::Operator::OperatorParameters>
    getParameters( const std::string &type,
                   AMP::LinearAlgebra::Vector::const_shared_ptr u,
                   AMP::shared_ptr<AMP::Operator::OperatorParameters> params = nullptr ) override;

protected:
    LinearTimeOperator();

    void getFromInput( const AMP::shared_ptr<AMP::Database> &db );

    bool d_bModifyRhsOperatorMatrix;

    /**
     * set to true if this operator corresponds to an algebraic component
     */
    bool d_bAlgebraicComponent;

    double d_dScalingFactor;
    double d_dCurrentDt;

    AMP::shared_ptr<AMP::Operator::LinearOperator> d_pRhsOperator;

    AMP::shared_ptr<AMP::Operator::LinearOperator> d_pMassOperator;

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pPreviousTimeSolution;

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pScratchVector;

    double d_current_time;
    double d_beta;

private:
};
} // namespace TimeIntegrator
} // namespace AMP

#endif
