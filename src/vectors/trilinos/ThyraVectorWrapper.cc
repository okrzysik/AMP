#include "vectors/trilinos/ThyraVectorWrapper.h"

namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Constructors                                                  *
****************************************************************/
ThyraVectorWrapper::ThyraVectorWrapper( AMP::LinearAlgebra::Vector::shared_ptr vec )
{
    d_vec = vec;
}


/****************************************************************
* Destructor                                                    *
****************************************************************/
ThyraVectorWrapper::~ThyraVectorWrapper( )
{
}


/****************************************************************
* Public functions derived from Thyra::LinearOpBase             *
****************************************************************/
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ThyraVectorWrapper::range() const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<const Thyra::VectorSpaceBase<double> >();
}
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ThyraVectorWrapper::domain() const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<const Thyra::VectorSpaceBase<double> >();
}

/****************************************************************
* Public functions derived from Thyra::MultiVectorBase          *
****************************************************************/
Teuchos::RCP<Thyra::MultiVectorBase<double> > ThyraVectorWrapper::clone_mv() const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<Thyra::MultiVectorBase<double> >();
}

/****************************************************************
* Public functions derived from Thyra::VectorBase               *
****************************************************************/
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ThyraVectorWrapper::space() const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<const Thyra::VectorSpaceBase<double> >();
}
Teuchos::RCP<Thyra::VectorBase<double> > ThyraVectorWrapper::clone_v() const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<Thyra::VectorBase<double> >();
}

/****************************************************************
* Protected functions derived from Thyra::LinearOpBase          *
****************************************************************/
bool ThyraVectorWrapper::opSupportedImpl(Thyra::EOpTransp M_trans) const
{
    AMP_ERROR("Not finished");
    return false;
}
void ThyraVectorWrapper::applyImpl(const Thyra::EOpTransp M_trans, const Thyra::MultiVectorBase<double> &X, 
        const Teuchos::Ptr<Thyra::MultiVectorBase<double> > &Y, const double alpha, const double beta) const
{
    AMP_ERROR("Not finished");
}

/****************************************************************
* Protected functions derived from Thyra::MultiVectorBase       *
****************************************************************/
Teuchos::RCP<Thyra::VectorBase<double> > ThyraVectorWrapper::nonconstColImpl(Teuchos::Ordinal j)
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<Thyra::VectorBase<double> >();
}
Teuchos::RCP<const Thyra::MultiVectorBase<double> > ThyraVectorWrapper::contigSubViewImpl(const Teuchos::Range1D &colRng) const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<const Thyra::MultiVectorBase<double> >();
}
Teuchos::RCP<Thyra::MultiVectorBase<double> > ThyraVectorWrapper::nonconstContigSubViewImpl(const Teuchos::Range1D &colRng)
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<Thyra::MultiVectorBase<double> >();
}
Teuchos::RCP<const Thyra::MultiVectorBase<double> > ThyraVectorWrapper::nonContigSubViewImpl(const Teuchos::ArrayView< const int > &cols) const
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<const Thyra::MultiVectorBase<double> >();
}
Teuchos::RCP<Thyra::MultiVectorBase<double> > ThyraVectorWrapper::nonconstNonContigSubViewImpl(const Teuchos::ArrayView< const int > &cols)
{
    AMP_ERROR("Not finished");
    return Teuchos::RCP<Thyra::MultiVectorBase<double> >();
}
void ThyraVectorWrapper::mvMultiReductApplyOpImpl( const RTOpPack::RTOpT<double> &primary_op, 
        const Teuchos::ArrayView< const Teuchos::Ptr< const Thyra::MultiVectorBase<double> > > &multi_vecs, 
        const Teuchos::ArrayView< const Teuchos::Ptr< Thyra::MultiVectorBase<double> > > &targ_multi_vecs, 
        const Teuchos::ArrayView< const Teuchos::Ptr< RTOpPack::ReductTarget > > &reduct_objs, 
        const Teuchos::Ordinal primary_global_offset ) const
{
    AMP_ERROR("Not finished");
}
void ThyraVectorWrapper::mvSingleReductApplyOpImpl( const RTOpPack::RTOpT<double> &primary_op, 
        const RTOpPack::RTOpT<double> &secondary_op, 
        const Teuchos::ArrayView< const Teuchos::Ptr< const Thyra::MultiVectorBase<double> > > &multi_vecs, 
        const Teuchos::ArrayView< const Teuchos::Ptr< Thyra::MultiVectorBase<double> > > &targ_multi_vecs, 
        const Teuchos::Ptr< RTOpPack::ReductTarget > &reduct_obj, 
        const Teuchos::Ordinal primary_global_offset ) const
{
    AMP_ERROR("Not finished");
}
void ThyraVectorWrapper::acquireDetachedMultiVectorViewImpl( const Teuchos::Range1D &rowRng, const Teuchos::Range1D &colRng, 
        RTOpPack::ConstSubMultiVectorView<double> *sub_mv ) const
{
    AMP_ERROR("Not finished");
}
void ThyraVectorWrapper::releaseDetachedMultiVectorViewImpl(RTOpPack::ConstSubMultiVectorView<double> *sub_mv) const
{
    AMP_ERROR("Not finished");
}
void ThyraVectorWrapper::acquireNonconstDetachedMultiVectorViewImpl(const Teuchos::Range1D &rowRng, const Teuchos::Range1D &colRng, 
        RTOpPack::SubMultiVectorView<double> *sub_mv )
{
    AMP_ERROR("Not finished");
}
void ThyraVectorWrapper::commitNonconstDetachedMultiVectorViewImpl( RTOpPack::SubMultiVectorView<double> *sub_mv )
{
    AMP_ERROR("Not finished");
}

/****************************************************************
* Protected functions derived from Thyra::VectorBase            *
****************************************************************/
void ThyraVectorWrapper::applyOpImpl( const RTOpPack::RTOpT<double> &op, 
        const Teuchos::ArrayView<const Teuchos::Ptr<const Thyra::VectorBase<double> > > &vecs, 
        const Teuchos::ArrayView<const Teuchos::Ptr<Thyra::VectorBase<double> > > &targ_vecs, 
        const Teuchos::Ptr<RTOpPack::ReductTarget> &reduct_obj, 
        const Teuchos::Ordinal global_offset) const
{
    AMP_ERROR("Not finished");
}
void ThyraVectorWrapper::acquireDetachedVectorViewImpl(const Teuchos::Range1D &rng, RTOpPack::ConstSubVectorView<double> *sub_vec) const
{
    AMP_ERROR("Not finished");
}
void ThyraVectorWrapper::releaseDetachedVectorViewImpl(RTOpPack::ConstSubVectorView<double> *sub_vec) const
{
    AMP_ERROR("Not finished");
}
void ThyraVectorWrapper::acquireNonconstDetachedVectorViewImpl(const Teuchos::Range1D &rng, RTOpPack::SubVectorView<double> *sub_vec)
{
    AMP_ERROR("Not finished");
}
void ThyraVectorWrapper::commitNonconstDetachedVectorViewImpl(RTOpPack::SubVectorView<double> *sub_vec)
{
    AMP_ERROR("Not finished");
}
void ThyraVectorWrapper::setSubVectorImpl(const RTOpPack::SparseSubVectorT<double> &sub_vec)
{
    AMP_ERROR("Not finished");
}


}
}

