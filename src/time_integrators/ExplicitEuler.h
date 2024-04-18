//
// $Id: ExplicitEuler.h,v 1.2 2006/02/07 17:36:01 philipb Exp $
// $Revision: 1.2 $
// $Date: 2006/02/07 17:36:01 $
//
// File:  ExplicitEuler.h
// Copyright:  (c) 2005 The Regents of the University of California
// Description:  Concrete time integrator using backward Euler method
//

#ifndef included_ExplicitEuler
#define included_ExplicitEuler

#include <string>

namespace AMP::LinearAlgebra {
class Vector;
}

#include "AMP/time_integrators/TimeIntegrator.h"

namespace AMP::TimeIntegrator {

/** \class ExplicitEuler
 *
 * Class ExplicitEuler is a concrete time integrator
 * that implements the explicit Runge-Kutta second order (RK2) method.
 */
class ExplicitEuler : public AMP::TimeIntegrator::TimeIntegrator
{
public:
    /**
     * Constructor that accepts parameter list.
     */
    explicit ExplicitEuler(
        std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters );

    /**
     * Destructor.
     */
    ~ExplicitEuler();

    static std::unique_ptr<TimeIntegrator> createTimeIntegrator(
        std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters )
    {
        return std::unique_ptr<TimeIntegrator>( new ExplicitEuler( parameters ) );
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

    std::string type() const override { return "ExplicitEuler"; }

private:
    /**
     * Constructor.
     */
    ExplicitEuler() = delete;

    /**
     * setup the vectors used by BE
     */
    void setupVectors( void );

    std::shared_ptr<AMP::LinearAlgebra::Vector> d_new_solution = nullptr;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_f_vec        = nullptr;
};
} // namespace AMP::TimeIntegrator

#endif
