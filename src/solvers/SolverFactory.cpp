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
#include "solvers/SolverFactory.h"

#include "solvers/SolverStrategy.h"
#include "solvers/SolverStrategyParameters.h"

#ifdef USE_EXT_HYPRE
#include "solvers/hypre/BoomerAMGSolver.h"
#endif

#ifdef USE_TRILINOS_ML
#include "solvers/trilinos/ml/TrilinosMLSolver.h"
#endif

#ifdef USE_TRILINOS_MUELU
#include "solvers/trilinos/muelu/TrilinosMueLuSolver.h"
#endif

#ifdef USE_EXT_PETSC
#include "solvers/petsc/PetscKrylovSolver.h"
#endif

#include "solvers/BiCGSTABSolver.h"
#include "solvers/CGSolver.h"
#include "solvers/GMRESSolver.h"
#include "solvers/QMRCGSTABSolver.h"
#include "solvers/TFQMRSolver.h"

namespace AMP {
namespace Solver {

// register all known solver factories
void registerSolverFactories()
{
    auto &solverFactory = SolverFactory::getFactory();

#ifdef USE_TRILINOS_MUELU
    solverFactory.registerFactory( "TrilinosMueLuSolver", TrilinosMueLuSolver::createSolver );
#endif

#ifdef USE_TRILINOS_ML
    solverFactory.registerFactory( "TrilinosMLSolver", TrilinosMLSolver::createSolver );
#endif

#ifdef USE_EXT_HYPRE
    solverFactory.registerFactory( "BoomerAMGSolver", BoomerAMGSolver::createSolver );
#endif

#ifdef USE_EXT_PETSC
    solverFactory.registerFactory( "PetscKrylovSolver", PetscKrylovSolver::createSolver );
#endif

    solverFactory.registerFactory( "CGSolver", CGSolver::createSolver );
    solverFactory.registerFactory( "GMRESSolver", GMRESSolver::createSolver );
    solverFactory.registerFactory( "BiCGSTABSolver", BiCGSTABSolver::createSolver );
    solverFactory.registerFactory( "TFQMRSolver", TFQMRSolver::createSolver );
    solverFactory.registerFactory( "QMRCGSTABSolver", QMRCGSTABSolver::createSolver );
}
} // namespace Solver
} // end namespace AMP
