#include "vectors/trilinos/ThyraVectorWrapper.h"
#include "vectors/trilinos/ThyraVectorSpaceWrapper.h"


#include "RTOpPack_Types.hpp"
#include "RTOpPack_RTOpT_decl.hpp"
#include "RTOpPack_SPMD_apply_op_def.hpp"
#ifdef USE_EXT_MPI
    #include "Teuchos_DefaultMpiComm.hpp"
#else
    #include "Teuchos_DefaultSerialComm.hpp"
#endif


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Constructors                                                  *
****************************************************************/
ThyraVectorWrapper::ThyraVectorWrapper( AMP::LinearAlgebra::Vector::shared_ptr vec )
{
    d_vec = vec;
    #ifdef USE_EXT_MPI
        Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > mpi_wrapper( 
            new Teuchos::OpaqueWrapper<MPI_Comm>(d_vec->getComm().getCommunicator()) );
        d_comm = new Teuchos::MpiComm<RTOpPack::index_type>( mpi_wrapper );
    #else
        d_comm = new Teuchos::SerialComm<RTOpPack::index_type>();
    #endif

}


/****************************************************************
* Destructor                                                    *
****************************************************************/
ThyraVectorWrapper::~ThyraVectorWrapper( )
{
    delete d_comm;
}


/****************************************************************
* Public functions derived from Thyra::LinearOpBase             *
****************************************************************/
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ThyraVectorWrapper::range() const
{
    return space();
}
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ThyraVectorWrapper::domain() const
{
    return space();
}


/****************************************************************
* Public functions derived from Thyra::MultiVectorBase          *
****************************************************************/
Teuchos::RCP<Thyra::MultiVectorBase<double> > ThyraVectorWrapper::clone_mv() const
{
    return clone_v();
}


/****************************************************************
* Public functions derived from Thyra::VectorBase               *
****************************************************************/
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ThyraVectorWrapper::space() const
{
    return Teuchos::RCP<const Thyra::VectorSpaceBase<double> >( 
        new ThyraVectorSpaceWrapper( d_vec->getDOFManager() ) );
}
Teuchos::RCP<Thyra::VectorBase<double> > ThyraVectorWrapper::clone_v() const
{
    return Teuchos::RCP<Thyra::VectorBase<double> >( new ThyraVectorWrapper( d_vec->cloneVector() ) );
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
void ThyraVectorWrapper::acquireDetachedVectorViewImpl(const Teuchos::Range1D &rng, RTOpPack::ConstSubVectorView<double> *sub_vec) const
{
    Vector::const_shared_ptr vec = getVec();
    AMP_ASSERT( vec->getLocalSize()>0 );
    size_t lower = rng.lbound();
    size_t upper = std::min<size_t>(rng.ubound(),vec->getLocalSize());
    AMP_ASSERT( upper>lower );
    size_t size = upper-lower;
    double* ptr = new double[size];
    Teuchos::ArrayRCP<double> array(ptr,rng.lbound(),rng.lbound()+size,true);
    size_t *indices = new size_t[size];
    for (size_t i=0; i<size; i++)
        indices[i] = lower+i;
    vec->getValuesByLocalID(size,indices,ptr);
    delete [] indices;
    sub_vec = new RTOpPack::ConstSubVectorView<double>( array );
}
void ThyraVectorWrapper::releaseDetachedVectorViewImpl(RTOpPack::ConstSubVectorView<double> *sub_vec) const
{
}
void ThyraVectorWrapper::acquireNonconstDetachedVectorViewImpl(const Teuchos::Range1D &rng, RTOpPack::SubVectorView<double> *sub_vec)
{
    Vector::shared_ptr vec = getVec();
    AMP_ASSERT( vec->getLocalSize()>0 );
    size_t lower = rng.lbound();
    size_t upper = std::min<size_t>(rng.ubound(),vec->getLocalSize());
    AMP_ASSERT( upper>lower );
    size_t size = upper-lower;
    double* ptr = new double[size];
    Teuchos::ArrayRCP<double> array(ptr,rng.lbound(),rng.lbound()+size,true);
    size_t *indices = new size_t[size];
    for (size_t i=0; i<size; i++)
        indices[i] = lower+i;
    vec->getValuesByLocalID(size,indices,ptr);
    delete [] indices;
    sub_vec = new RTOpPack::SubVectorView<double>( array );
}
void ThyraVectorWrapper::commitNonconstDetachedVectorViewImpl(RTOpPack::SubVectorView<double> *sub_vec)
{
    Vector::shared_ptr vec = getVec();
    const Teuchos::ArrayRCP<double> &array = sub_vec->values();
    size_t lower = array.lowerOffset();
    size_t upper = array.upperOffset();
    AMP_ASSERT( upper>lower && upper<=vec->getLocalSize() );
    size_t size = upper-lower;
    size_t *indices = new size_t[size];
    for (size_t i=0; i<size; i++)
        indices[i] = lower+i;
    vec->setValuesByLocalID(size,indices,array.get());
    delete [] indices;
}
void ThyraVectorWrapper::setSubVectorImpl(const RTOpPack::SparseSubVectorT<double> &sub_vec)
{
    AMP_ERROR("Not finished");
}


/****************************************************************
* applyOpImpl derived from Thyra::VectorBase                    *
****************************************************************/
void ThyraVectorWrapper::applyOpImpl( const RTOpPack::RTOpT<double> &op, 
        const Teuchos::ArrayView<const Teuchos::Ptr<const Thyra::VectorBase<double> > > &vecs, 
        const Teuchos::ArrayView<const Teuchos::Ptr<Thyra::VectorBase<double> > > &targ_vecs, 
        const Teuchos::Ptr<RTOpPack::ReductTarget> &reduct_obj, 
        const Teuchos::Ordinal global_offset) const
{
    size_t n_blocks = d_vec->numberOfDataBlocks();
    std::vector<size_t> block_size(n_blocks,0);
    for (size_t i=0; i<n_blocks; i++)
        block_size[i] = d_vec->sizeOfDataBlock(i);
    // Check that all vectors are compatible
    for (int i=0; i<vecs.size(); i++) {
        const ThyraVectorWrapper* ptr = dynamic_cast<const ThyraVectorWrapper*>( vecs[i].getRawPtr() );
        AMP_INSIST(ptr!=NULL,"All vectors used in applyOpImpl must be of the type ThyraVectorWrapper");
        AMP_ASSERT(ptr->d_vec->numberOfDataBlocks()==n_blocks);
        for (size_t j=0; j<n_blocks; j++)
            AMP_ASSERT(block_size[j]==ptr->d_vec->sizeOfDataBlock());
    }
    for (int i=0; i<targ_vecs.size(); i++) {
        const ThyraVectorWrapper* ptr = dynamic_cast<const ThyraVectorWrapper*>( targ_vecs[i].getRawPtr() );
        AMP_INSIST(ptr!=NULL,"All vectors used in applyOpImpl must be of the type ThyraVectorWrapper");
        AMP_ASSERT(ptr->d_vec->numberOfDataBlocks()==n_blocks);
        for (size_t j=0; j<n_blocks; j++)
            AMP_ASSERT(block_size[j]==ptr->d_vec->sizeOfDataBlock());
    }
    // Apply the operation
    Teuchos::RCP<RTOpPack::ReductTarget> reduct_obj2 = op.reduct_obj_create();
    for (size_t j=0; j<n_blocks; j++) {
        std::vector<RTOpPack::ConstSubVectorView<double> >  sub_vecs(vecs.size());
        std::vector<RTOpPack::SubVectorView<double> > targ_sub_vecs(targ_vecs.size());
        for (int i=0; i<vecs.size(); i++) {
            const ThyraVectorWrapper* ptr = dynamic_cast<const ThyraVectorWrapper*>( vecs[i].getRawPtr() );
            sub_vecs[i] = RTOpPack::ConstSubVectorView<double>( 
                Teuchos::ArrayRCP<double>(ptr->d_vec->getRawDataBlock<double>(j),0,block_size[j],false) );
        }
        for (int i=0; i<targ_vecs.size(); i++) {
            const ThyraVectorWrapper* ptr = dynamic_cast<const ThyraVectorWrapper*>( targ_vecs[i].getRawPtr() );
            targ_sub_vecs[i] = RTOpPack::SubVectorView<double>( 
                Teuchos::ArrayRCP<double>(ptr->d_vec->getRawDataBlock<double>(j),0,block_size[j],false) );
        }
        Teuchos::ArrayView<const RTOpPack::ConstSubVectorView<double> > sub_vecs2 = 
             Teuchos::ArrayView<const RTOpPack::ConstSubVectorView<double> >(sub_vecs);
        Teuchos::ArrayView<const RTOpPack::SubVectorView<double> > targ_sub_vecs2 = 
           Teuchos::ArrayView<const RTOpPack::SubVectorView<double> >(targ_sub_vecs);
        op.apply_op( sub_vecs2, targ_sub_vecs2, reduct_obj2.ptr() );
    }
    // Reduce the result
    if ( reduct_obj.get()!=NULL ) {
        RTOpPack::SPMD_all_reduce<double>( d_comm, op, 1, 
            Teuchos::tuple<const RTOpPack::ReductTarget*>(&*reduct_obj2).getRawPtr(),
            Teuchos::tuple<RTOpPack::ReductTarget*>(&*reduct_obj).getRawPtr() );
    }
    // Change the vector state targ_vecs
    for (int i=0; i<targ_vecs.size(); i++) {
        ThyraVectorWrapper* ptr = dynamic_cast<ThyraVectorWrapper*>( targ_vecs[i].get() );
        if ( ptr!=NULL ) {
            if ( ptr->d_vec->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED )
                ptr->d_vec->setUpdateStatus( AMP::LinearAlgebra::Vector::LOCAL_CHANGED );
        }
    }
}


}
}

