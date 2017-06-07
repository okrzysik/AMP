
#ifndef included_AMP_TrilinosMueLuSolver
#define included_AMP_TrilinosMueLuSolver

#include "matrices/trilinos/EpetraMatrix.h"
#include "solvers/SolverStrategy.h"
#include "solvers/SolverStrategyParameters.h"
#include "solvers/trilinos/ml/MLoptions.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include "Teuchos_ParameterList.hpp"
// Xpetra include
#include <Xpetra_Parameters.hpp>
#include <Xpetra_Operator_fwd.hpp>
#include "MueLu_FactoryManager_decl.hpp"
//#include "MuelLu_Hierarchy.hpp"
//#include <MueLu.hpp>
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop

namespace MueLu{
  class EpetraOperator;
  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node> class Hierarchy;
  
  using Scalar=double;
  using LocalOrdinal=int;
  using GlobalOrdinal=int;
  using Node=Xpetra::EpetraNode;
}

namespace AMP {
namespace Solver {


using TrilinosMueLuSolverParameters = SolverStrategyParameters;
using SC=MueLu::Scalar;
using LO=MueLu::LocalOrdinal;
using GO=MueLu::GlobalOrdinal;
using NO=MueLu::Node;

/**
 * The TrilinosMueLuSolver is a wrapper to the Trilinos ML solver. ML provides implementations of
 * various algebraic multigrid methods. The wrapper at present simply provides an adaptor
 * to enable AMP users to use the black box ML preconditioner.
 */

class TrilinosMueLuSolver : public SolverStrategy
{

public:
    /**
     * Default constructor
     */
    TrilinosMueLuSolver();

    /**
     * Main constructor.
     @param [in] parameters The parameters object contains a database object which must contain the
     following fields in addition to the fields expected by the base class SolverStrategy class:

    */
    explicit TrilinosMueLuSolver( AMP::shared_ptr<TrilinosMueLuSolverParameters> parameters );

    /**
     * Default destructor
     */
    virtual ~TrilinosMueLuSolver();

    /**
     * Solve the system \f$Au = f\f$.
     @param [in] f : shared pointer to right hand side vector
     @param [out] u : shared pointer to approximate computed solution
     */
    void solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                AMP::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

    /**
     * Solve the system \f$Au = f\f$.
     @param [in] f : shared pointer to right hand side vector
     @param [out] u : shared pointer to approximate computed solution
     */
    void solveWithHierarchy( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                             AMP::shared_ptr<AMP::LinearAlgebra::Vector> u );

    /**
     * Return a shared pointer to the ML_Epetra::MultiLevelPreconditioner object
     */
    inline const AMP::shared_ptr<MueLu::EpetraOperator> getMLSolver( void )
    {
        return d_mueluSolver;
    }

    /**
     * Initialize the solution vector and potentially create internal vectors needed for solution
     @param [in] parameters The parameters object
     contains a database object. Refer to the documentation for the constructor to see what fields
     are required.
     This routine assumes that a non-NULL operator of type LinearOperator has been registered with
     the solver.
     The LinearOperator currently is assumed to contain a pointer to an EpetraMatrix object.
     */
    void initialize( AMP::shared_ptr<SolverStrategyParameters> const parameters ) override;

    /**
     * Register the operator that the solver will use during solves
     @param [in] op shared pointer to the linear operator $A$ for equation \f$A u = f\f$
     */
    void registerOperator( const AMP::shared_ptr<AMP::Operator::Operator> op ) override;

    /**
     * Resets the associated operator internally with new parameters if necessary
     * @param [in] params
     *        OperatorParameters object that is NULL by default
     */
    void resetOperator( const AMP::shared_ptr<AMP::Operator::OperatorParameters> params ) override;

    /**
     * Resets the solver internally with new parameters if necessary
     * @param [in] params
     *        SolverStrategyParameters object that is NULL by default
     * Currently every call to reset destroys the ML preconditioner object
     * and recreates it based on the parameters object. See constructor for
     * fields required for parameter object.
     */
    void reset( AMP::shared_ptr<SolverStrategyParameters> params ) override;

protected:
    void reSolveWithLU( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                        AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                            u );

    void getFromInput( const AMP::shared_ptr<AMP::Database> &db );

private:
    bool d_bUseEpetra;
    bool d_build_from_components = false;  //! whether to explicitly build the hierarchy
    AMP_MPI d_comm;

    bool d_bCreationPhase; /**< set to true if the PC is not ready and false otherwise. */
    bool d_bRobustMode;

    AMP::shared_ptr<MueLu::EpetraOperator> d_mueluSolver;

    AMP::shared_ptr<AMP::LinearAlgebra::EpetraMatrix> d_matrix;
    Teuchos::ParameterList d_MueLuParameterList;

    Teuchos::RCP< MueLu::Hierarchy<SC,LO,GO,NO> > d_mueluHierarchy; //! AMG hierarchy

    MueLu::FactoryManager<SC,LO,GO,NO> d_factoryManager; //! factory manager for MueLu components

};
}
}

#endif
