#ifndef included_AMP_ExplicitEuler
#define included_AMP_ExplicitEuler

#include <string>

#ifndef included_AMP_TimeIntegrator
#include "TimeIntegrator.h"
#endif

namespace AMP {
namespace TimeIntegrator {

/** \class ExplicitEuler
 *
 * Class ExplicitEuler is a concrete time integrator
 * that implements the explicit Runge-Kutta second order (RK2) method.
 */
class ExplicitEuler : public TimeIntegrator
{
public:
    /**
     * Constructor that accepts parameter list.
     */
    explicit ExplicitEuler( AMP::shared_ptr<TimeIntegratorParameters> parameters );

    /**
     * Destructor.
     */
    virtual ~ExplicitEuler();

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
    ExplicitEuler();

    /**
     * Read data from input database.
     */
    void getFromInput( AMP::shared_ptr<AMP::Database> input_db );

    /**
     * setup the vectors used by BE
     */
    void setupVectors( void );

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_new_solution;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_f_vec;
};
} // namespace TimeIntegrator
} // namespace AMP

#endif
