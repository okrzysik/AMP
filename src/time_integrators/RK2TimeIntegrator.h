//
// $Id: RK2TimeIntegrator.h,v 1.2 2006/02/07 17:36:01 philipb Exp $
// $Revision: 1.2 $
// $Date: 2006/02/07 17:36:01 $
//
// File:  RK2TimeIntegrator.h
// Copyright:  (c) 2005 The Regents of the University of California
// Description:  Concrete time integrator using backward Euler method
//

#ifndef included_RK2TimeIntegrator
#define included_RK2TimeIntegrator

#include <string>

#include "AMP/time_integrators/TimeIntegrator.h"

namespace AMP::TimeIntegrator {

class TimeIntegratorParameters;

/** \class RK2TimeIntegrator
 *
 * Class RK2TimeIntegrator is a concrete time integrator
 * that implements the explicit Runge-Kutta second order (RK2) method
 * also known as Heun's method.
 */
class RK2TimeIntegrator : public AMP::TimeIntegrator::TimeIntegrator
{
public:
    /**
     * Constructor that accepts parameter list.
     */
    explicit RK2TimeIntegrator( std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> );

    /**
     * Destructor.
     */
    ~RK2TimeIntegrator();

    static std::unique_ptr<AMP::TimeIntegrator::TimeIntegrator> createTimeIntegrator(
        std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters )
    {
        return std::unique_ptr<AMP::TimeIntegrator::TimeIntegrator>(
            new RK2TimeIntegrator( parameters ) );
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

private:
    /**
     * Constructor.
     */
    RK2TimeIntegrator() = delete;

    /**
     * setup the vectors used by RK2
     */
    void setupVectors( void );

    int d_number_regrid_states = 0;

    std::shared_ptr<AMP::LinearAlgebra::Vector> d_new_solution = nullptr;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_k1_vec       = nullptr;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_k2_vec       = nullptr;
};
} // namespace AMP::TimeIntegrator

#endif
