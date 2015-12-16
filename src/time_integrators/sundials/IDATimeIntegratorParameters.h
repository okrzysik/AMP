#ifndef included_IDATimeIntegratorParameters
#define included_IDATimeIntegratorParameters


#ifndef included_TimeIntegratorParameters
#include "time_integrators/TimeIntegratorParameters.h"
#endif

#ifndef included_IDATimeOperator
#include "time_integrators/sundials/IDATimeOperator.h"
#endif

#ifndef included_LinearTimeOperator
#include "time_integrators/LinearTimeOperator.h"
#endif

// BC : Altered this to get a compile going..
//#ifndef included_AMP_MassMatrix
//#include "operators/MassMatrix.h"
//#endif
#include "operators/Operator.h"
#include "operators/libmesh/MassLinearFEOperator.h"

#ifndef included_SolverStrategy
#include "solvers/SolverStrategy.h"
#endif

namespace AMP {
namespace TimeIntegrator {

typedef AMP::Operator::MassLinearFEOperator MassOperator;

/*!
 @brief TimeIntegratorParameters is a base class for providing
 parameters for the TimeIntegrator's. The Database object contained
 must contain the following entries:

 Required input keys and data types:
 @param initial_time double value for the initial simulation time.
 @param final_time double value for the final simulation time.
 @param max_integrator_steps integer value for the maximum number
 of timesteps allowed.

 All input data items described above, except for initial_time,
 may be overwritten by new input values when continuing from restart.

 */

class IDATimeIntegratorParameters : public TimeIntegratorParameters {
public:
    explicit IDATimeIntegratorParameters( const AMP::shared_ptr<AMP::Database> db );

    virtual ~IDATimeIntegratorParameters();


    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_ic_vector_prime;

    // Needs to be fixed - JL
    AMP::shared_ptr<IDATimeOperator> d_pIDATimeOperator;
    AMP::shared_ptr<LinearTimeOperator> d_pLinearTimeOperator;
    AMP::shared_ptr<TimeOperatorParameters> d_pLinearTimeOperatorParameters;

    AMP::shared_ptr<AMP::Solver::SolverStrategy> d_pPreconditioner;
    AMP::shared_ptr<AMP::Operator::LinearOperator> d_pLinearOperator;


    // source term
    // temp operators
    AMP::shared_ptr<AMP::Operator::Operator> d_temp_operator_1;
    AMP::shared_ptr<AMP::Operator::Operator> d_temp_operator_2;

protected:
private:
    // not implemented
    // IDATimeIntegratorParameters(){}

    explicit IDATimeIntegratorParameters(); // Just following ImplicitTimeIntegratorParameters();
    explicit IDATimeIntegratorParameters( const IDATimeIntegratorParameters & );
    void operator=( const IDATimeIntegratorParameters & );
};
}
}

#endif
