#ifndef included_AMP_contact_MPCSolver
#define included_AMP_contact_MPCSolver

#include "AMP/solvers/SolverStrategy.h"

namespace AMP::Solver {

typedef SolverStrategyParameters MPCSolverParameters;

class MPCSolver : public SolverStrategy
{
public:
    explicit MPCSolver( std::shared_ptr<MPCSolverParameters> params ) : SolverStrategy( params ) {}

    std::string type() const override { return "MPCSolver"; }

    void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                std::shared_ptr<AMP::LinearAlgebra::Vector> u );
};
} // namespace AMP::Solver

#endif
