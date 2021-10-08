//
// $Id: RK45TimeIntegrator.h,v 1.2 2006/02/07 17:36:01 philipb Exp $
// $Revision: 1.2 $
// $Date: 2006/02/07 17:36:01 $
//
// File:  RK45TimeIntegrator.h
// Copyright:  (c) 2005 The Regents of the University of California
// Description:  Concrete time integrator using Runge-Kutta RK4 method
//

#ifndef included_RK45TimeIntegrator
#define included_RK45TimeIntegrator

#include <string>

#include "AMP/time_integrators/TimeIntegrator.h"

namespace AMP::TimeIntegrator {

/** \class RK45TimeIntegrator
 *
 * Class RK45TimeIntegrator is a concrete time integrator
 * that implements the explicit Runge-Kutta fourth order (RK4) method.
 */
class RK45TimeIntegrator : public AMP::TimeIntegrator::TimeIntegrator
{
public:
    /**
     * Constructor that accepts parameter list.
     */
    explicit RK45TimeIntegrator(
        std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters );

    /**
     * Destructor.
     */
    ~RK45TimeIntegrator();

    static std::unique_ptr<AMP::TimeIntegrator::TimeIntegrator> createTimeIntegrator(
        std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters )
    {
        return std::unique_ptr<AMP::TimeIntegrator::TimeIntegrator>(
            new RK45TimeIntegrator( parameters ) );
    }

    /**
     * Initialize from parameter list.
     */
    virtual void initialize(
        std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters ) override;

    /**
     * Resets the internal state of the time integrator as needed.
     * A parameter argument is passed to allow for general flexibility
     * in determining what needs to be reset Typically used after a regrid.
     */
    virtual void reset(
        std::shared_ptr<const AMP::TimeIntegrator::TimeIntegratorParameters> parameters ) override;

    /**
     * Specify next time step to use.
     */
    virtual double getNextDt( const bool good_solution ) override;

    /**
     * Determine whether time advanced solution is satisfactory.
     */
    virtual bool checkNewSolution( void ) override;

    /**
     * Update state of the solution.
     */
    virtual void updateSolution( void ) override;

    virtual int advanceSolution( const double dt,
                                 const bool first_step,
                                 std::shared_ptr<AMP::LinearAlgebra::Vector> in,
                                 std::shared_ptr<AMP::LinearAlgebra::Vector> out ) override;

private:
    /**
     * Constructor.
     */
    RK45TimeIntegrator();

    /**
     * Read data from input database.
     */
    void getFromInput( std::shared_ptr<AMP::Database> input_db );

    /**
     * setup the vectors used by RK4
     */
    void setupVectors( void );

    int d_number_regrid_states = 0;

    double d_safety_factor = 0.0;
    double d_atol          = 0.0;

    bool d_use_fixed_dt = false;

    std::shared_ptr<AMP::LinearAlgebra::Vector> d_new_solution = nullptr;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_k1_vec       = nullptr;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_k2_vec       = nullptr;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_k3_vec       = nullptr;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_k4_vec       = nullptr;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_k5_vec       = nullptr;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_k6_vec       = nullptr;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_z_vec        = nullptr;
};
} // namespace AMP::TimeIntegrator

#endif
