#include "solvers/CGSolver.h"
#include "vectors/petsc/ManagedPetscVector.h"
#include "operators/LinearOperator.h"
#include "matrices/Matrix.h"
#include "matrices/petsc/PetscMatrix.h"
#include "ProfilerApp.h"

extern "C"{
#include "assert.h"
}

namespace AMP {
namespace Solver {

/****************************************************************
*  Constructors                                                 *
****************************************************************/
CGSolver::CGSolver()
{
}
CGSolver::CGSolver(AMP::shared_ptr<CGSolverParameters> parameters):SolverStrategy(parameters)
{
    assert(parameters.get()!=NULL);

    // Initialize
    initialize(parameters);
}


/****************************************************************
*  Destructor                                                   *
****************************************************************/
CGSolver::~CGSolver()
{
}

/****************************************************************
*  Initialize                                                   *
****************************************************************/
void
CGSolver::initialize(AMP::shared_ptr<SolverStrategyParameters> const params)
{
    AMP::shared_ptr<CGSolverParameters> parameters = AMP::dynamic_pointer_cast<CGSolverParameters>(params);
    AMP_ASSERT(parameters.get()!=NULL);
    d_comm = parameters->d_comm;
    AMP_ASSERT(!d_comm.isNull());

    d_pPreconditioner = parameters->d_pPreconditioner;

    getFromInput(parameters->d_db);

    if(d_pOperator.get()!=NULL) {
      registerOperator(d_pOperator);
    }
}

// Function to get values from input
void CGSolver::getFromInput(const AMP::shared_ptr<AMP::Database> &db)
{

  d_dRelativeTolerance = db->getDoubleWithDefault("relative_tolerance", 1.0e-9);
  d_dAbsoluteTolerance = db->getDoubleWithDefault("absolute_tolerance", 1.0e-14);
  d_dDivergenceTolerance = db->getDoubleWithDefault("divergence_tolerance", 1.0e+03);
  d_iMaxIterations       = db->getDoubleWithDefault("max_iterations", 1000);

  d_bUsesPreconditioner = db->getBoolWithDefault("uses_preconditioner", false);
}

/****************************************************************
*  Solve                                                        *
****************************************************************/
void
CGSolver::solve(AMP::shared_ptr<const AMP::LinearAlgebra::Vector>  f,
                  AMP::shared_ptr<AMP::LinearAlgebra::Vector>  u)
{
    PROFILE_START("solve");

    // Check input vector states
    AMP_ASSERT( (f->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED) ||
		(f->getUpdateStatus() == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
    AMP_ASSERT( (u->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED) ||
		(u->getUpdateStatus() == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
    
    if(d_iDebugPrintInfoLevel>1) {
      std::cout << "CGSolver::solve: initial L2Norm of solution vector: " << u->L2Norm() << std::endl;
      std::cout << "CGSolver::solve: initial L2Norm of rhs vector: " << f->L2Norm() << std::endl;
    }

    if(d_bUsesPreconditioner) {
    } else {
    }
    if(d_pOperator.get()!=NULL) {
       registerOperator(d_pOperator);
    }

    if(d_iDebugPrintInfoLevel>2) {
        std::cout << "L2Norm of solution: " << u->L2Norm() << std::endl;
    }

    PROFILE_STOP("solve");
}

/****************************************************************
*  Function to set the register the operator                    *
****************************************************************/
void CGSolver::registerOperator(const AMP::shared_ptr<AMP::Operator::Operator> op)
{
  // in this case we make the assumption we can access a PetscMat for now
  assert(op.get()!=NULL);

  d_pOperator = op;

  AMP::shared_ptr<AMP::Operator::LinearOperator> linearOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>(op);
  assert(linearOperator.get() != NULL);

}
void CGSolver::resetOperator(const AMP::shared_ptr<AMP::Operator::OperatorParameters> params)
{
  if(d_pOperator.get()!=NULL) {
    d_pOperator->reset(params);    
  }
  
  // should add a mechanism for the linear operator to provide updated parameters for the preconditioner operator
  // though it's unclear where this might be necessary
  if(d_pPreconditioner.get()!=NULL) {
    d_pPreconditioner->resetOperator(params);
  }
}

}
}

