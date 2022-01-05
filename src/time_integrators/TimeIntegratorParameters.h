#ifndef included_AMP_TimeIntegratorParameters
#define included_AMP_TimeIntegratorParameters

#include "AMP/operators/Operator.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/ParameterBase.h"
#include "AMP/vectors/Vector.h"
#include <memory>

#include <string>


namespace AMP::TimeIntegrator {

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

class TimeIntegratorParameters : public Operator::OperatorParameters
{
public:
    //! Convience typedef
    typedef std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> shared_ptr;

    explicit TimeIntegratorParameters( std::shared_ptr<AMP::Database> db );

    virtual ~TimeIntegratorParameters();

    /**
     * String used to identify specific class instantiation of the time integrator
     */
    std::string d_object_name;

    AMP::AMP_MPI d_comm; // Comm for this object

    /**
     * Initial conditions vector
     */
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_ic_vector;

    /**
     * source term for time integration, can also include boundary conditions for IBVP problems
     */
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_pSourceTerm;

    /**
     * The operator is the right hand side operator for an explicit integrator when the time
     * integration problem is :
     * u_t = f(u)
     * but in the case of implicit time integrators the operator represents u_t-f(u)
     */
    std::shared_ptr<AMP::Operator::Operator> d_operator;

    /**
     * The operator is the left hand side mass operator (for FEM formulations)
     */
    std::shared_ptr<AMP::Operator::Operator> d_pMassOperator;

    /**
     * algebraic variable
     */
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_pAlgebraicVariable;

    //! pointer to global database
    std::shared_ptr<AMP::Database> d_global_db;

protected:
private:
    // not implemented
    TimeIntegratorParameters()                                            = delete;
    explicit TimeIntegratorParameters( const TimeIntegratorParameters & ) = delete;
    void operator=( const TimeIntegratorParameters & ) = delete;
};
} // namespace AMP::TimeIntegrator

#endif
