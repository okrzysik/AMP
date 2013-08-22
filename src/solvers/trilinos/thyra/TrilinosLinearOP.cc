#include "solvers/trilinos/thyra/TrilinosLinearOP.h"
#include "vectors/trilinos/thyra/ThyraVector.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Solver {


/****************************************************************
*  Constructors                                                 *
****************************************************************/
TrilinosLinearOP::TrilinosLinearOP()
{
}
TrilinosLinearOP::TrilinosLinearOP( AMP::Operator::Operator::shared_ptr op )
{
    this->d_operator = op;
    AMP_ASSERT(d_operator!=NULL);
}
TrilinosLinearOP::TrilinosLinearOP( AMP::Solver::SolverStrategy::shared_ptr solver )
{
    this->d_solver = solver;
    AMP_ASSERT(d_solver!=NULL);
}
TrilinosLinearOP::~TrilinosLinearOP()
{
}


/****************************************************************
*  Functions inherited from Thyra::LinearOpBase                 *
****************************************************************/
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > TrilinosLinearOP::range() const
{
    AMP_ERROR("Not Implimented Yet");
    return Teuchos::RCP<const Thyra::VectorSpaceBase<double> >();
}
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > TrilinosLinearOP::domain () const
{
    AMP_ERROR("Not Implimented Yet");
    return Teuchos::RCP<const Thyra::VectorSpaceBase<double> >();
}
bool TrilinosLinearOP::opSupportedImpl(Thyra::EOpTransp) const
{
    AMP_ERROR("Not Implimented Yet");
    return false;
}
void TrilinosLinearOP::applyImpl(const Thyra::EOpTransp M_trans, const Thyra::MultiVectorBase<double> &X, 
        const Teuchos::Ptr< Thyra::MultiVectorBase<double> > &Y, const double alpha, const double beta) const
{
    // Compute Y = alpha*OP(M)*X + beta*Y 
    AMP_ASSERT(M_trans==Thyra::NOTRANS);
    AMP::LinearAlgebra::Vector::const_shared_ptr x = 
        AMP::LinearAlgebra::ThyraVector::constView( dynamic_cast<const Thyra::VectorBase<double>*>(&X) );
    AMP::LinearAlgebra::Vector::shared_ptr y = 
        AMP::LinearAlgebra::ThyraVector::view( dynamic_cast<Thyra::VectorBase<double>*>(Y.get()) );
    //x->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    //y->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    if ( d_operator != NULL ) {
        // Apply the AMP::Operator to compute the residual
        AMP::LinearAlgebra::Vector::shared_ptr f = y->cloneVector();
        d_operator->apply( f, x, y, alpha, beta );
    } else {
        // Apply the AMP::Solver
        d_solver->solve( x, y );
    }
}


}
}

