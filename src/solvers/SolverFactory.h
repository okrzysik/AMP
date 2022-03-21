#ifndef included_AMP_SolverFactory_H
#define included_AMP_SolverFactory_H

#include "AMP/utils/FactoryStrategy.hpp"

namespace AMP::Solver {

class SolverStrategy;
class SolverStrategyParameters;

using SolverFactory = FactoryStrategy<SolverStrategy, SolverStrategyParameters>;

// free function to preregister solvers known by AMP
void registerSolverFactories();
} // namespace AMP::Solver
#endif
