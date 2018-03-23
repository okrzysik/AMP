#ifndef included_AMP_RK4TimeIntegrator
#define included_AMP_RK4TimeIntegrator

#include <string>

#ifndef included_AMP_TimeIntegrator
#include "TimeIntegrator.h"
#endif

namespace AMP {
namespace TimeIntegrator {

/** \class RK4TimeIntegrator
 *
 * Class RK4TimeIntegrator is a concrete time integrator
 * that implements the explicit Runge-Kutta fourth order (RK4) method.
 */
class RK4TimeIntegrator : public TimeIntegrator
{
public:
    /**
     * Constructor that accepts parameter list.
     */
    explicit RK4TimeIntegrator( AMP::shared_ptr<TimeIntegratorParameters> parameters );

    /**
     * Destructor.
     */
    virtual ~RK4TimeIntegrator();

    /**
     * Initialize from parameter list.
     */
    void initialize( AMP::shared_ptr<TimeIntegratorParameters> parameters ) override;

    /**
     * Resets the internal state of the time integrator as needed.
     * A parameter argument is passed to allow for general flexibility
     * in determining what needs to be reset Typically used after a regrid.
     */
    void reset( AMP::shared_ptr<TimeIntegratorParameters> parameters ) override;

    /**
     * Specify initial time step.
     */
    double getInitialDt();

    /**
     * Specify next time step to use.
     */
    double getNextDt( const bool good_solution ) override;

    /**
     * Determine whether time advanced solution is satisfactory.
     */
    bool checkNewSolution( void ) const override;

    /**
     * Update state of the solution.
     */
    void updateSolution( void ) override;

    int advanceSolution( const double dt, const bool first_step ) override;

private:
    /**
     * Constructor.
     */
    RK4TimeIntegrator();

    /**
     * Read data from input database.
     */
    void getFromInput( AMP::shared_ptr<AMP::Database> input_db );

    /**
     * setup the vectors used by RK4
     */
    void setupVectors( void );

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_new_solution;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_k1_vec;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_k2_vec;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_k3_vec;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_k4_vec;
};
} // namespace TimeIntegrator
} // namespace AMP

#endif
