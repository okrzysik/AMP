#include "solvers/TrilinosLinearOP.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Solver {


/****************************************************************
*  Constructors                                                 *
****************************************************************/
TrilinosLinearOP::TrilinosLinearOP()
{
    
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
    AMP_ERROR("Not Implimented Yet");
}


}
}

