#ifndef included_AMP_FlowMapOperator
#define included_AMP_FlowMapOperator

#include "operators/map/Map3to1to3.h"

namespace AMP {
namespace Operator {

class FlowMapOperator : public Map3to1to3 {
public:
    // Overload applyStart
    virtual void applyStart( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                             AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr r,
                             const double a = -1.0,
                             const double b = 1.0 );

    // Overload applyFinish
    virtual void applyFinish( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                              AMP::LinearAlgebra::Vector::const_shared_ptr u,
                              AMP::LinearAlgebra::Vector::shared_ptr r,
                              const double a = -1.0,
                              const double b = 1.0 );

protected:
    // Internal data
    std::multimap<double, double> d_TempClad, d_TempFlow;
    double d_K, d_De, d_Re, d_Pr, d_G, Cp;

    // Overload buildMap
    virtual void buildMap( const AMP::LinearAlgebra::Vector::shared_ptr p );
    virtual void buildMap( const AMP::LinearAlgebra::Vector::shared_ptr p,
                           std::multimap<double, double> &m );

    // Overload buildReturn
    virtual void buildReturn( AMP::LinearAlgebra::Vector::shared_ptr p );

    // Overload smear
    virtual void smear( double tolerance );
    virtual void smear( double tolerance, std::multimap<double, double> &m );
};
}
}

#endif
