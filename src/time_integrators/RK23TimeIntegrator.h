//
// $Id: RK23TimeIntegrator.h,v 1.2 2006/02/07 17:36:01 philipb Exp $
// $Revision: 1.2 $
// $Date: 2006/02/07 17:36:01 $
//
// File:  RK23TimeIntegrator.h
// Copyright:  (c) 2005 The Regents of the University of California
// Description:  Concrete time integrator using Runge-Kutta RK23 method
//

#ifndef included_RK23TimeIntegrator
#define included_RK23TimeIntegrator

#include <string>

#include "AMP/time_integrators/TimeIntegrator.h"

namespace AMP::TimeIntegrator {

/** \class RK23TimeIntegrator
 *
 * Class RK23TimeIntegrator is a concrete time integrator
 * that implements the explicit Bogacki-Shampine adaptive Runge-Kutta (Matlab
 * ode23) method.
 */
class RK23TimeIntegrator : public AMP::TimeIntegrator::TimeIntegrator
{
public:
    /**
     * Constructor that accepts parameter list.
     */
    explicit RK23TimeIntegrator(
        std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters );

    /**
     * Destructor.
     */
    ~RK23TimeIntegrator();

    static std::unique_ptr<AMP::TimeIntegrator::TimeIntegrator> createTimeIntegrator(
        std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters )
    {
        return std::unique_ptr<AMP::TimeIntegrator::TimeIntegrator>(
            new RK23TimeIntegrator( parameters ) );
    }

    /**
     * Initialize from parameter list.
     */
    void initialize(
        std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters ) override;

    /**
     * Resets the internal state of the time integrator as needed.
     * A parameter argument is passed to allow for general flexibility
     * in determining what needs to be reset Typically used after a regrid.
     */
    void reset(
        std::shared_ptr<const AMP::TimeIntegrator::TimeIntegratorParameters> parameters ) override;

    /**
     * Specify next time step to use.
     */
    double getNextDt( const bool good_solution ) override;

    /**
     * Determine whether time advanced solution is satisfactory.
     */
    bool checkNewSolution( void ) override;

    /**
     * Update state of the solution.
     */
    void updateSolution( void ) override;

    int advanceSolution( const double dt,
                         const bool first_step,
                         std::shared_ptr<AMP::LinearAlgebra::Vector> in,
                         std::shared_ptr<AMP::LinearAlgebra::Vector> out ) override;

    std::string type() const override { return "RK23"; }

private:
    /**
     * Constructor.
     */
    RK23TimeIntegrator();

    /**
     * Read data from input database.
     */
    void getFromInput( std::shared_ptr<AMP::Database> input_db );

    /**
     * setup the vectors used by RK23
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
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_z_vec        = nullptr;
};
} // namespace AMP::TimeIntegrator

#endif
