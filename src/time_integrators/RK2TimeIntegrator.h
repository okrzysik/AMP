#ifndef included_AMP_RK2TimeIntegrator
#define included_AMP_RK2TimeIntegrator

#include <string>

#ifndef included_AMP_TimeIntegrator
#include "TimeIntegrator.h"
#endif

namespace AMP {
namespace TimeIntegrator {

/** \class RK2TimeIntegrator
 *
 * Class RK2TimeIntegrator is a concrete time integrator
 * that implements the explicit Runge-Kutta second order (RK2) method.
 */
class RK2TimeIntegrator : public TimeIntegrator
{
public:
    /**
     * Constructor that accepts parameter list.
     */
    explicit RK2TimeIntegrator( std::shared_ptr<TimeIntegratorParameters> parameters );

    /**
     * Destructor.
     */
    virtual ~RK2TimeIntegrator();

    /**
     * Initialize from parameter list.
     */
    void initialize( std::shared_ptr<TimeIntegratorParameters> parameters ) override;

    /**
     * Resets the internal state of the time integrator as needed.
     * A parameter argument is passed to allow for general flexibility
     * in determining what needs to be reset Typically used after a regrid.
     */
    void reset( std::shared_ptr<TimeIntegratorParameters> parameters ) override;

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
    RK2TimeIntegrator();

    /**
     * Read data from input database.
     */
    void getFromInput( std::shared_ptr<AMP::Database> input_db );

    /**
     * setup the vectors used by RK2
     */
    void setupVectors( void );

    std::shared_ptr<AMP::LinearAlgebra::Vector> d_new_solution;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_k1_vec;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_k2_vec;
};
} // namespace TimeIntegrator
} // namespace AMP

#endif
