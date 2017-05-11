#ifndef included_IDATimeOperator
#define included_IDATimeOperator


#include "operators/Operator.h"
#include "operators/OperatorBuilder.h"
#include "operators/OperatorParameters.h"
#include "time_integrators/TimeOperator.h"
#include "time_integrators/TimeOperatorParameters.h"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"
#include "vectors/Vector.h"


namespace AMP {
namespace TimeIntegrator {

typedef TimeOperatorParameters IDATimeOperatorParameters;

/*!
  @brief operator class associated with IDATimeIntegrator

  Class IDATimeOperator is derived from TimeOperator. It
  is the operator class associated with a IDATimeIntegrator.

  @see IDATimeIntegrator
  @see TimeOperator
*/

class IDATimeOperator : public TimeOperator
{
public:
    /**
     * Main constructor.
     @param [in] params: shared pointer to TimeOperatorParameters object.
     */
    explicit IDATimeOperator( AMP::shared_ptr<AMP::Operator::OperatorParameters> params );

    /**
     * virtual destructor
     */
    virtual ~IDATimeOperator();

    // virtual void reset(const AMP::shared_ptr<AMP::Operator::OperatorParameters>& params);

    /**
      The function that computes the residual.
     * @param u: multivector of the state.
     * @param f: The result of apply ( f = A(u) )
     */
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    /**
     * registers the time derivative vector provided by IDA with this operator
     @param [in] vec   shared pointer to time derivative computed by IDA
     */
    void registerIDATimeDerivative( AMP::shared_ptr<AMP::LinearAlgebra::Vector> vec )
    {
        d_pIDATimeDerivative = vec;
    }

    /**
     * registers a source term if any
     @param [in] vec   shared pointer to vector for source term
     */
    void registerSourceTerm( AMP::shared_ptr<AMP::LinearAlgebra::Vector> vec )
    {
        d_pSourceTerm = vec;
    }

    /**
     * sets the current time for the operator
     @param [in] currentTime   the current time
     */
    void registerCurrentTime( double currentTime ) { d_current_time = currentTime; }

protected:
    IDATimeOperator();

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pIDATimeDerivative;

    bool d_cloningHappened;

    // JL
    // The test we want to run has a source term which depends on time
    // The time comes from TimeIntegrator
    double d_current_time;
    double d_beta;

private:
};
}
}

#endif
