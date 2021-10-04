#ifndef included_TimeIntegratorInterface_h_
#define included_TimeIntegratorInterface_h_

#include <memory>
#include <string>
#include <vector>

namespace AMP::LinearAlgebra {
class Vector;
}

namespace AMP::TimeIntegrator {
/**
 * This class provides the interface required by the BDF methods at present
 */
class TimeIntegratorInterface
{
public:
    TimeIntegratorInterface()          = default;
    virtual ~TimeIntegratorInterface() = default;

    /**
     * Compute r=f(u) where u_t = f(u)
     *  \pre It is assumed that all coarse fine ghost values have been set
     * appropriately
     * \param u
     *           vector for solution variables
     * \param r
     *           vector for residual
     */
    virtual void applyRhs( std::shared_ptr<const AMP::LinearAlgebra::Vector> x,
                           std::shared_ptr<AMP::LinearAlgebra::Vector> f ) = 0;

    virtual void setTimeOperatorScaling( const double gamma ) { d_gamma = gamma; }

    virtual double getGamma( void ) const { return d_gamma; }

protected:
    //! stores linear weighted sum of previous time solutions that becomes a source term
    //! for the nonlinear system of equations
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_integrator_source_vector = nullptr;

    double d_gamma = 0.0;
};
} // namespace AMP::TimeIntegrator

#endif
