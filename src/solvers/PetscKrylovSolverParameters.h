#ifndef included_AMP_PetscKrylovSolverParameters
#define included_AMP_PetscKrylovSolverParameters

#ifndef included_Pointer
#include "boost/shared_ptr.hpp"
#endif

#ifndef included_Database
#include "utils/Database.h"
#endif

#ifndef included_AMP_SolverStrategyParameters
#include "SolverStrategyParameters.h"
#endif

extern "C"{
#ifdef MPICH_SKIP_MPICXX
#define _FIX_FOR_PETSC_MPI_CXX
#undef MPICH_SKIP_MPICXX
#endif

#ifdef OMPI_SKIP_MPICXX
#define _FIX_FOR_PETSC_OMPI_CXX
#undef OMPI_SKIP_MPICXX
#endif

#include "petsc.h"

#ifdef _FIX_FOR_PETSC_OMPI_CXX
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif
#endif

#ifdef _FIX_FOR_PETSC_MPI_CXX
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#endif

}

#ifndef included_AMP_SolverStrategy
#include "SolverStrategy.h"
#endif

namespace AMP {
namespace Solver {

  /**
   * Class PetscKrylovSolverParameters provides a uniform mechanism to pass
   * initialization parameters to the PetscKrylovSolver solver. It contains
   * shared pointers to a database object and a preconditioner. All member variables are public.
   */
  class PetscKrylovSolverParameters: public SolverStrategyParameters{
  public:
    PetscKrylovSolverParameters(){}
    PetscKrylovSolverParameters(const boost::shared_ptr<AMP::Database> db);
    ~PetscKrylovSolverParameters(){}

    AMP_MPI d_comm;

    boost::shared_ptr<AMP::Solver::SolverStrategy> d_pPreconditioner;

  protected:
  private:
    
  };

}
}

#endif
