#include "solvers/GMRESSolver.h"
#include "operators/LinearOperator.h"
#include "ProfilerApp.h"

extern "C"{
#include "assert.h"
}

#include <cmath>
#include <limits>

namespace AMP {
namespace Solver {

/****************************************************************
*  Constructors                                                 *
****************************************************************/
GMRESSolver::GMRESSolver(): d_restarts(0)
{
}

GMRESSolver::GMRESSolver(AMP::shared_ptr<KrylovSolverParameters> parameters):SolverStrategy(parameters), d_restarts(0)
{
    assert(parameters.get()!=NULL);
    
    // Initialize
    initialize(parameters);
}


/****************************************************************
*  Destructor                                                   *
****************************************************************/
GMRESSolver::~GMRESSolver()
{
}

/****************************************************************
*  Initialize                                                   *
****************************************************************/
void
GMRESSolver::initialize(AMP::shared_ptr<SolverStrategyParameters> const params)
{
    AMP::shared_ptr<KrylovSolverParameters> parameters = AMP::dynamic_pointer_cast<KrylovSolverParameters>(params);
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
void GMRESSolver::getFromInput(const AMP::shared_ptr<AMP::Database> &db)
{

  d_dRelativeTolerance = db->getDoubleWithDefault("relative_tolerance", 1.0e-9);
  d_iMaxIterations       = db->getDoubleWithDefault("max_iterations", 1000);

  d_bUsesPreconditioner = db->getBoolWithDefault("uses_preconditioner", false);
}

/****************************************************************
*  Solve                                                        *
* TODO: store convergence history, iterations, convergence reason
****************************************************************/
void
GMRESSolver::solve(AMP::shared_ptr<const AMP::LinearAlgebra::Vector>  f,
		   AMP::shared_ptr<AMP::LinearAlgebra::Vector>  u)
{
  PROFILE_START("solve");
  
  // Check input vector states
  AMP_ASSERT( (f->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED) ||
	      (f->getUpdateStatus() == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
  AMP_ASSERT( (u->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED) ||
	      (u->getUpdateStatus() == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );

  // compute the norm of the rhs in order to compute
  // the termination criterion
  double f_norm = f->L2Norm();

  // if the rhs is zero we try to converge to the relative convergence
  if(f_norm == 0.0 ) {
    f_norm = 1.0;
  }

  const double terminate_tol = d_dRelativeTolerance*f_norm;
  
  if(d_iDebugPrintInfoLevel>1) {
    std::cout << "GMRESSolver::solve: initial L2Norm of solution vector: " << u->L2Norm() << std::endl;
    std::cout << "GMRESSolver::solve: initial L2Norm of rhs vector: " << f_norm << std::endl;
  }

  if(d_pOperator.get()!=NULL) {
    registerOperator(d_pOperator);
  }
     
  // residual vector
  AMP::LinearAlgebra::Vector::shared_ptr res = f->cloneVector();
  
  // compute the initial residual
  if( d_bUseZeroInitialGuess ) {
    res->copyVector(f);
  } else {
    d_pOperator->residual(f, u, res);
  }

  // compute the current residual norm
  double res_norm = res->L2Norm();

  // exit if the residual is already low enough
  if(res_norm < terminate_tol ) {
    // provide a convergence reason
    // provide history (iterations, conv history etc)
    return;
  }

  if(d_iDebugPrintInfoLevel>2) {
    std::cout << "L2Norm of solution: " << u->L2Norm() << std::endl;
  }
  
  PROFILE_STOP("solve");
}

/****************************************************************
*  Function to set the register the operator                    *
****************************************************************/
void GMRESSolver::registerOperator(const AMP::shared_ptr<AMP::Operator::Operator> op)
{
  assert(op.get()!=NULL);

  d_pOperator = op;

  AMP::shared_ptr<AMP::Operator::LinearOperator> linearOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>(op);
  assert(linearOperator.get() != NULL);

}
void GMRESSolver::resetOperator(const AMP::shared_ptr<AMP::Operator::OperatorParameters> params)
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

