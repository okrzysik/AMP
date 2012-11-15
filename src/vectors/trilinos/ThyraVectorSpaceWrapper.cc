#include "vectors/trilinos/ThyraVectorSpaceWrapper.h"

#include "utils/Utilities.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Constructors                                                  *
****************************************************************/
ThyraVectorSpaceWrapper::ThyraVectorSpaceWrapper( AMP::Discretization::DOFManager::shared_ptr dofs )
{
    d_dofs = dofs;
}


/****************************************************************
* Destructor                                                    *
****************************************************************/
ThyraVectorSpaceWrapper::~ThyraVectorSpaceWrapper( )
{
}


/****************************************************************
* Virtual functions inherited from VectorSpaceBase              *
****************************************************************/
Teuchos::Ordinal ThyraVectorSpaceWrapper::dim() const
{
    return d_dofs->numGlobalDOF();
}
bool ThyraVectorSpaceWrapper::isCompatible(const Thyra::VectorSpaceBase<double> &vecSpc) const
{
    AMP_ERROR("Not finished");
    return false;
}
Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<double> > ThyraVectorSpaceWrapper::smallVecSpcFcty() const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<double> >();
}
double ThyraVectorSpaceWrapper::scalarProd(const Thyra::VectorBase<double> &x, const Thyra::VectorBase<double> &y) const
{
    AMP_ERROR("Not finished");
    return 0.0;
}
Teuchos::RCP<Thyra::VectorBase<double> >  ThyraVectorSpaceWrapper::createMember() const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<Thyra::VectorBase<double> >();
}
Teuchos::RCP<Thyra::MultiVectorBase <double> >  ThyraVectorSpaceWrapper::createMembers(int numMembers) const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<Thyra::MultiVectorBase<double> >();
}
Teuchos::RCP<Thyra::VectorBase<double> >  ThyraVectorSpaceWrapper::createMemberView(const RTOpPack::SubVectorView<double> &raw_v) const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<Thyra::VectorBase<double> >();
}
Teuchos::RCP<const Thyra::VectorBase<double> >  ThyraVectorSpaceWrapper::createMemberView(const RTOpPack::ConstSubVectorView<double> &raw_v) const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<const Thyra::VectorBase<double> >();
}
Teuchos::RCP<Thyra::MultiVectorBase<double> >  ThyraVectorSpaceWrapper::createMembersView(const RTOpPack::SubMultiVectorView<double> &raw_mv) const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<Thyra::MultiVectorBase<double> >();
}
Teuchos::RCP<const Thyra::MultiVectorBase<double> >  ThyraVectorSpaceWrapper::createMembersView(const RTOpPack::ConstSubMultiVectorView<double> &raw_mv) const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<const Thyra::MultiVectorBase<double> >();
}
void ThyraVectorSpaceWrapper::scalarProdsImpl(const Thyra::MultiVectorBase<double> &X, const Thyra::MultiVectorBase<double> &Y, const Teuchos::ArrayView<double> &scalarProds) const
{
    AMP_ERROR("Not finished");
}


}
}

