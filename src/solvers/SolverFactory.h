/*
Copyright 2005, The Regents of the University
of California. This software was produced under
a U.S. Government contract (W-7405-ENG-36)
by Los Alamos National Laboratory, which is
operated by the University of California for the
U.S. Department of Energy. The U.S.
Government is licensed to use, reproduce, and
distribute this software. Permission is granted
to the public to copy and use this software
without charge, provided that this Notice and
any statement of authorship are reproduced on
all copies. Neither the Government nor the
University makes any warranty, express or
implied, or assumes any liability or
responsibility for the use of this software.
*/

#ifndef included_SolverFactory_H_
#define included_SolverFactory_H_

#include "utils/FactoryStrategy.hpp"

namespace AMP {

namespace Solver {

class SolverStrategy;
class SolverStrategyParameters;

using SolverFactory = FactoryStrategy<SolverStrategy, SolverStrategyParameters>;

// free function to preregister solvers known by SAMRSolvers
void registerSolverFactories();
} // namespace Solver
} // namespace AMP
#endif // included_SolverFactory_H_
