#ifndef included_TimeOperator
#define included_TimeOperator

#ifndef included_AMP_config

#endif

#ifndef included_AMP_Operator
#include "operators/Operator.h"
#endif

#ifndef included_Pointer
#include "utils/shared_ptr.h"
#endif

#ifndef included_AMP_Vector
#include "vectors/Vector.h"
#endif

#include "utils/Utilities.h"

#ifndef included_AMP_OperatorParameters
#include "operators/OperatorParameters.h"
#endif

// BP the following two lines are probably unnecessary, not removing
// prior to Jan, 2011 code camp.
// JL
#include "operators/libmesh/VolumeIntegralOperator.h"

namespace AMP {
namespace TimeIntegrator {

/*!
  @brief base class for operator class associated with ImplicitTimeIntegrator

  Class ImplicitTimeOperator is a base class derived from Operator. It
  is the operator class associated with a ImplicitTimeIntegrator. The solver associated
  with the ImplicitTimeIntegrator will register this object.
  For implicit time integrators solving a problem M(u_t) = f(u), at each time step
  a nonlinear or linear problem, F(u) = 0, where F(u) = M(u_t)-f(u) discretized, needs to be solved.
  The TimeOperator class is meant to represent this operator. Internally, two operators are used to
  represent
  F, a mass operator and a rhs operator. M corresponds to the mass operator and f() corresponds to
  the rhs operator.
  The mass operator is NULL in the case of FVM, FD discretizations.

  @see ImplicitTimeIntegrator
  @see Operator
  @see SolverStrategy
  @see TimeOperatorParameters
*/
class TimeOperator : public AMP::Operator::Operator
{
public:
    /**
     * Main constructor meant for use by users.
     @param [in] params : parameter object that must be of type
     TimeOperatorParameters. The database object contained in the parameter
     object contains the following fields:

     1. name: bLinearMassOperator
        description: boolean to indicate whether the mass operator is a linear operator, currently
        either both mass and rhs operators have to be linear or both have to be nonlinear
        type: bool
        default value: FALSE
        valid values: (TRUE, FALSE)
        optional field: yes

     2. name: bLinearRhsOperator
        description: boolean to indicate whether the rhs operator is a linear operator, currently
        either both mass and rhs operators have to be linear or both have to be nonlinear
        type: bool
        default value: FALSE
        valid values: (TRUE, FALSE)
        optional field: yes

     3. name: bAlgebraicComponent
        description: for a DAE system this TimeOperator would be one of several components. This
     field
        indicates if this component is an algebraic component, ie, having no time derivative
        type: bool
        default value: FALSE
        valid values: (TRUE, FALSE)
        optional field: yes

     4. name: CurrentDt
        description: current time step value
        type: double
        default value: none
        valid values: (positve real values)
        optional field: no

    */
    explicit TimeOperator( AMP::shared_ptr<AMP::Operator::OperatorParameters> params );

    /**
     * virtual destructor
     */
    virtual ~TimeOperator();

    /**
     * This function is useful for re-initializing an operator
     * \param params
     *        parameter object containing parameters to change
     */
    virtual void reset(const AMP::shared_ptr<AMP::Operator::OperatorParameters>
                           &params) override;

    /**
     * This function registers a rhs operator with the TimeOperator class
     @param [in] op : shared pointer to Operator, cannot be another TimeOperator
     */
    void registerRhsOperator( AMP::shared_ptr<AMP::Operator::Operator> op ) { d_pRhsOperator = op; }

    /**
     * This function registers a mass operator with the TimeOperator class. Not necessary
     * for FD or FVM discretizations
     @param [in] op : shared pointer to Operator, cannot be another TimeOperator
     */
    void registerMassOperator( AMP::shared_ptr<AMP::Operator::Operator> op )
    {
        d_pMassOperator = op;
    }

    /**
     * register a variable as being an algebraic component. Deprecated.
     */
    void registerAlgebraicVariable( AMP::shared_ptr<AMP::LinearAlgebra::Variable> var )
    {
        d_pAlgebraicVariable = var;
    }

    /**
     * return a shared pointer to the rhs operator
     */
    AMP::shared_ptr<AMP::Operator::Operator> getRhsOperator( void ) { return d_pRhsOperator; }

    /**
     * return a shared pointer to the mass operator
     */
    AMP::shared_ptr<AMP::Operator::Operator> getMassOperator( void ) { return d_pMassOperator; }

    /**
     * register a vector consisting of the solution at the previous time step.
     @param [in] previousSolution : shared pointer to Vector
     */
    void setPreviousSolution( AMP::shared_ptr<AMP::LinearAlgebra::Vector> previousSolution )
    {
        d_pPreviousTimeSolution = previousSolution;
    }

    /**
     * set the value of the current time step
     @param [in] dt : value of current timestep
     */
    void setDt( double dt ) { d_dCurrentDt = dt; }

    /**
     * returns a Variable object corresponding to the output variable for the TimeOperator.
     * Currently
     * this is the output variable associated with the rhs operator.
     */
    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override {
      return d_pRhsOperator->getOutputVariable();
    }

    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    /**
     * implements the getJacobianParameters interface required by operators. This routine returns
     * a shared pointer to a TimeOperatorParameters object, internally containing shares pointerst
     to
     * two OperatorParameter objects, one for the mass operator (if not a FD or FVM discretization)
     * and one for the rhs operator.
     @param [in] u : shared pointer to a Vector at which the Jacobian is to be evaluated.
     */
    AMP::shared_ptr<AMP::Operator::OperatorParameters>
    getParameters( const std::string &type,
                   AMP::LinearAlgebra::Vector::const_shared_ptr u,
                   AMP::shared_ptr<AMP::Operator::OperatorParameters> params = NULL ) override;

protected:
    TimeOperator();

    void getFromInput( const AMP::shared_ptr<AMP::Database> &db );

    bool d_bLinearMassOperator;

    bool d_bLinearRhsOperator;

    /**
     * set to true if this operator corresponds to an algebraic component
     */
    bool d_bAlgebraicComponent;

    double d_dCurrentDt;

    /**
     * pointer to rhs operator
     */
    AMP::shared_ptr<AMP::Operator::Operator> d_pRhsOperator;

    /**
     * pointer to mass operator
     */
    AMP::shared_ptr<AMP::Operator::Operator> d_pMassOperator;

    /**
     * algebraic variable
     */
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_pAlgebraicVariable;

    /**
     * solution at previous time step
     */
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pPreviousTimeSolution;

    /**
     * scratch vector for internal use
     */
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pScratchVector;

    /**
     * vector containing source terms if any
     */
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pSourceTerm;

private:
};
}
}

#endif
