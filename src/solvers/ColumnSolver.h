#ifndef included_ColumnSolver
#define included_ColumnSolver

#include <string>
#include <vector>
#include "SolverStrategy.h"

namespace AMP {
namespace Solver {

  typedef SolverStrategyParameters ColumnSolverParameters;

  /**
   * A class for representing a column of solvers to be applied to a ColumnOperator
   * The main intention of this class is to provide block diagonal or block triangular
   * preconditioners of Gauss-Seidel or Jacobi type
   */
  class ColumnSolver : public SolverStrategy
  {
    public:

    /**
     * default empty constructor
     */
      ColumnSolver(){}

      /**
       * Main constructor
       @param [in] parameters object used for initialization. The parameters object
       contains a database object. The database object can optionally contain a
       string field "IterationType" in addition to any fields expected by the base
       class SolverStrategy. The IterationType must currently either be "GaussSeidel"
       or "SymmetricGaussSeidel"
       */
      ColumnSolver(boost::shared_ptr<SolverStrategyParameters> parameters);

      /**
       * destructor, currently does nothing
       */
      virtual ~ColumnSolver(){}

      /**
       * Solve the system \f$A(u) = f\f$.
       Assumptions: $A(.)$ is assumed to be a ColumnOperator
       @param [in] f shared pointer to right hand side
       @param [out] u shared pointer to computed approximate solution
       */
      void solve(boost::shared_ptr<AMP::LinearAlgebra::Vector>  f, boost::shared_ptr<AMP::LinearAlgebra::Vector>  u);

      /**
       * sets the initial guess
       @param [in] initialGuess shared pointer to initialGuess vector
       */
      void setInitialGuess( boost::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess );

      /**
       * @param [in] solver
       *            shared pointer to a solver to append to the existing column of solvers
       */
      void append(boost::shared_ptr<AMP::Solver::SolverStrategy> solver);

      /**
       * returns a shared pointer to the i-th solver
       @param [in] i integer index for solver to extract
       */
      boost::shared_ptr<AMP::Solver::SolverStrategy> getSolver(const int i){return d_Solvers[i]; }

      int getNumberOfSolvers() {
        return d_Solvers.size();
      }

      /**
       * Resets the associated operator internally with new parameters if necessary
       * @param parameters
       *        OperatorParameters object that is NULL by default
       */
      void resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> params);

    protected :

      void GaussSeidel(boost::shared_ptr<AMP::LinearAlgebra::Vector>  &f, boost::shared_ptr<AMP::LinearAlgebra::Vector>  &u);

      void SymmetricGaussSeidel(boost::shared_ptr<AMP::LinearAlgebra::Vector>  &f, boost::shared_ptr<AMP::LinearAlgebra::Vector>  &u);

      std::vector< boost::shared_ptr<AMP::Solver::SolverStrategy> > d_Solvers;
      /**
       * type of block iteration, valid values: GaussSeidel, SymmetricGaussSeidel
       */
      std::string d_IterationType;

      bool d_resetColumnOperator;
    private:
  };

}
}

#endif
