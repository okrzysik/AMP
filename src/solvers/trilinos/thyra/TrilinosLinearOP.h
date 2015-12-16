#ifndef included_AMP_TrilinosLinearOP
#define included_AMP_TrilinosLinearOP

// AMP includes
#include "operators/LinearOperator.h"
#include "solvers/SolverStrategy.h"
#include "vectors/Vector.h"

// Trilinos includes
#include "Thyra_LinearOpBase_def.hpp"


namespace AMP {
namespace Solver {


/**
  * The TrilinosLinearOP is a wrapper for a Thyra LinearOpBase to
  * wrap AMP::Operators or AMP::Solvers for use with Trilinos NOX solvers.
  */
class TrilinosLinearOP : public Thyra::LinearOpBase<double>
{
public:
    //! Constructor that wraps and AMP::Operator
    explicit TrilinosLinearOP( AMP::Operator::Operator::shared_ptr );

    //! Constructor that wraps and AMP::Solver
    explicit TrilinosLinearOP( AMP::Solver::SolverStrategy::shared_ptr );

    //! Destructor
    virtual ~TrilinosLinearOP();

    // Functions inherited from Thyra::LinearOpBase
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<double>> range() const;
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<double>> domain() const;
    virtual bool opSupportedImpl( Thyra::EOpTransp ) const;
    virtual void applyImpl( const Thyra::EOpTransp M_trans,
                            const Thyra::MultiVectorBase<double> &X,
                            const Teuchos::Ptr<Thyra::MultiVectorBase<double>> &Y,
                            const double alpha,
                            const double beta ) const;

private:
    //! Empty constructor
    TrilinosLinearOP();

    //! Data variables
    AMP::shared_ptr<AMP::Operator::Operator> d_operator;
    AMP::shared_ptr<AMP::Solver::SolverStrategy> d_solver;
};
}
}

#endif
