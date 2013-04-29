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
    d_op = boost::dynamic_pointer_cast<AMP::Operator::LinearOperator>( op );

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
    // Compute r = alpha*OP(M)*X + beta*Y 
    AMP_ASSERT(M_trans==Thyra::NOTRANS);
    AMP::LinearAlgebra::Vector::const_shared_ptr x = 
        AMP::LinearAlgebra::ThyraVector::constView( dynamic_cast<const Thyra::VectorBase<double>*>(&X) );
    AMP::LinearAlgebra::Vector::shared_ptr y = 
        AMP::LinearAlgebra::ThyraVector::view( dynamic_cast<Thyra::VectorBase<double>*>(Y.get()) );
    //x->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    //y->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    AMP::LinearAlgebra::Vector::shared_ptr f = y->cloneVector();
    d_op->apply( f, x, y, alpha, beta );
}


}
}

