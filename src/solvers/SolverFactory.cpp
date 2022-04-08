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


#include "AMP/solvers/SolverFactory.h"
#include "AMP/AMP_TPLs.h"
#include "AMP/solvers/BiCGSTABSolver.h"
#include "AMP/solvers/CGSolver.h"
#include "AMP/solvers/GMRESSolver.h"
#include "AMP/solvers/NonlinearKrylovAccelerator.h"
#include "AMP/solvers/QMRCGSTABSolver.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/TFQMRSolver.h"

#ifdef AMP_USE_PETSC
    #include "AMP/solvers/petsc/PetscKrylovSolver.h"
    #include "AMP/solvers/petsc/PetscSNESSolver.h"
#endif

#ifdef AMP_USE_HYPRE
    #include "AMP/solvers/hypre/BoomerAMGSolver.h"
#endif

#ifdef AMP_USE_TRILINOS_ML
    #include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#endif

#ifdef AMP_USE_TRILINOS_MUELU
    #include "AMP/solvers/trilinos/muelu/TrilinosMueLuSolver.h"
#endif


namespace AMP::Solver {

// register all known solver factories
void registerSolverFactories()
{
    auto &solverFactory = SolverFactory::getFactory();

#ifdef AMP_USE_TRILINOS_MUELU
    solverFactory.registerFactory( "TrilinosMueLuSolver", TrilinosMueLuSolver::createSolver );
#endif

#ifdef AMP_USE_TRILINOS_ML
    solverFactory.registerFactory( "TrilinosMLSolver", TrilinosMLSolver::createSolver );
#endif

#ifdef AMP_USE_HYPRE
    solverFactory.registerFactory( "BoomerAMGSolver", BoomerAMGSolver::createSolver );
#endif

#ifdef AMP_USE_PETSC
    solverFactory.registerFactory( "SNESSolver", PetscSNESSolver::createSolver );
    solverFactory.registerFactory( "PetscSNESSolver", PetscSNESSolver::createSolver );
    solverFactory.registerFactory( "PetscKrylovSolver", PetscKrylovSolver::createSolver );
#endif

    solverFactory.registerFactory( "CGSolver", CGSolver::createSolver );
    solverFactory.registerFactory( "GMRESSolver", GMRESSolver::createSolver );
    solverFactory.registerFactory( "BiCGSTABSolver", BiCGSTABSolver::createSolver );
    solverFactory.registerFactory( "TFQMRSolver", TFQMRSolver::createSolver );
    solverFactory.registerFactory( "QMRCGSTABSolver", QMRCGSTABSolver::createSolver );

    solverFactory.registerFactory( "NKASolver", NonlinearKrylovAccelerator::createSolver );
}


} // namespace AMP::Solver
