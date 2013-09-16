#include "vectors/trilinos/thyra/ThyraVectorWrapper.h"
#include "vectors/trilinos/thyra/ThyraVectorSpaceWrapper.h"


#include "RTOpPack_Types.hpp"
#include "RTOpPack_RTOpT_decl.hpp"
#include "RTOpPack_SPMD_apply_op_def.hpp"
#ifdef USE_EXT_MPI
    #include "Teuchos_DefaultMpiComm.hpp"
#else
    #include "Teuchos_DefaultSerialComm.hpp"
#endif
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"

namespace AMP {
namespace LinearAlgebra {


// Some function definitions
template<class TYPE>
static std::vector<ThyraVectorWrapper*> getPtr( std::vector<size_t> block_size,
    Teuchos::ArrayView<const Teuchos::Ptr<TYPE> > vecs );
template<class TYPE>
static std::vector<const ThyraVectorWrapper*> getConstPtr( std::vector<size_t> block_size,
    Teuchos::ArrayView<const Teuchos::Ptr<const TYPE> > vecs );



/****************************************************************
* Constructors                                                  *
****************************************************************/
ThyraVectorWrapper::ThyraVectorWrapper( const std::vector<AMP::LinearAlgebra::Vector::shared_ptr>& vecs )
{
    std::vector<size_t> cols(vecs.size());
    for (size_t i=0; i<vecs.size(); i++)
        cols[i] = i;
    this->initialize( vecs, cols, vecs.size() );
}
ThyraVectorWrapper::ThyraVectorWrapper( const std::vector<AMP::LinearAlgebra::Vector::shared_ptr>& vecs, 
    const std::vector<size_t>& cols, size_t N_cols )
{
    this->initialize( vecs, cols, N_cols );
}
void ThyraVectorWrapper::initialize( const std::vector<AMP::LinearAlgebra::Vector::shared_ptr>& vecs, 
    const std::vector<size_t>& cols, size_t N_cols )
{
    AMP_ASSERT(!vecs.empty());
    AMP_ASSERT(vecs.size()==cols.size());
    for (size_t i=0; i<vecs.size(); i++) {
        AMP_ASSERT(vecs[i]!=NULL);
        AMP_ASSERT(cols[i]<N_cols);
        for (size_t j=0; j<i; j++) {
            AMP_ASSERT(vecs[i]!=vecs[j]);
            AMP_ASSERT(cols[i]>cols[j]);
        }
    }
    // Check that the DOFs are compatible for all copies of the vector
    AMP::Discretization::DOFManager::shared_ptr dofs1 = vecs[0]->getDOFManager();
    for (size_t i=1; i<vecs.size(); i++) {
        AMP::Discretization::DOFManager::shared_ptr dofs2 = vecs[i]->getDOFManager();
        if ( dofs1==dofs2 )
            continue;
        if ( *dofs1!=*dofs2 )
            AMP_ERROR("The DOFManagers for all copies must match");
    }
    d_vecs = vecs;
    d_cols = cols;
    d_N_cols = N_cols;
    #ifdef USE_EXT_MPI
        Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > mpi_wrapper( 
            new Teuchos::OpaqueWrapper<MPI_Comm>(d_vecs[0]->getComm().getCommunicator()) );
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
* Get a shared object to this                                   *
****************************************************************/
template<class T> static void nullDeleter( T* ) {};
boost::shared_ptr<const ThyraVectorWrapper> ThyraVectorWrapper::shared_from_this() const
{
    return boost::shared_ptr<const ThyraVectorWrapper>( this, nullDeleter<const ThyraVectorWrapper> );
}


/****************************************************************
* Functions to return the number of rows and columns            *
****************************************************************/
size_t ThyraVectorWrapper::numRows() const
{
    // The number of rows is the number of unknowns in the vectors
    return d_vecs[0]->getDOFManager()->numGlobalDOF();
}
size_t ThyraVectorWrapper::numColumns() const
{
    // The number of columns is the number of internal vectors
    return d_N_cols;
}


/****************************************************************
* Functions to return the vector space                          *
****************************************************************/
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ThyraVectorWrapper::range() const
{
    boost::shared_ptr<const ThyraVectorWrapper> this_ptr = shared_from_this();
    return Teuchos::RCP<const Thyra::VectorSpaceBase<double> >( 
        new ThyraVectorSpaceWrapper( this_ptr, true ) );
}
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ThyraVectorWrapper::domain() const
{
    boost::shared_ptr<const ThyraVectorWrapper> this_ptr = shared_from_this();
    return Teuchos::RCP<const Thyra::VectorSpaceBase<double> >( 
        new ThyraVectorSpaceWrapper( this_ptr, false ) );
}
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ThyraVectorWrapper::space() const
{
    return range();
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
Teuchos::RCP<Thyra::VectorBase<double> > ThyraVectorWrapper::clone_v() const
{
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vecs(d_vecs.size());
    for (size_t i=0; i<d_vecs.size(); i++)
        vecs[i] = d_vecs[i]->cloneVector();
    return Teuchos::RCP<Thyra::VectorBase<double> >( new ThyraVectorWrapper( vecs, d_cols, d_N_cols ) );
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
    size_t M_dim[2], X_dim[2], Y_dim[2];
    M_dim[0] = this->numRows();                             // Number of rows of M
    M_dim[1] = this->numColumns();                          // Number of columns of M
    X_dim[0] = X.range()->dim();                            // Number of rows of X
    X_dim[1] = X.domain()->dim();                           // Number of columns of X
    Y_dim[0] = Y->range()->dim();                           // Number of rows of Y
    Y_dim[1] = Y->domain()->dim();                          // Number of columns of Y
    if ( M_trans==Thyra::NOTRANS || M_trans==Thyra::CONJ ) {
        // We are computing: Y = alpha*M*X + beta*Y 
        AMP_ASSERT( M_dim[1]==X_dim[0] && M_dim[0]==Y_dim[0] && X_dim[1]==Y_dim[1] );
        if ( X_dim[1]==1 ) {
            // We are performing a series of axpby operations
            const ThyraVectorWrapper* y = dynamic_cast<const ThyraVectorWrapper*>(Y.get());
            AMP_ASSERT(y!=NULL);
            AMP_ASSERT(d_cols.size()==y->d_cols.size());
            for (size_t i=0; i<d_cols.size(); i++)
                AMP_ASSERT(d_cols[i]==y->d_cols[i]);
            for (size_t i=0; i<d_cols.size(); i++) {
                double x = Thyra::get_ele(*(X.col(0)),d_cols[i]);
                y->d_vecs[i]->axpby(alpha*x,beta,d_vecs[i]);
            }
        } else {
            AMP_ERROR("Not finished");
        }
    } else if ( M_trans==Thyra::TRANS || M_trans==Thyra::CONJTRANS ) {
        // We are computing: Y = alpha*transpose(M)*X + beta*Y 
        AMP_ASSERT( M_dim[0]==X_dim[0] && M_dim[1]==Y_dim[0] && X_dim[1]==Y_dim[1] );
        if ( dynamic_cast<const ThyraVectorWrapper*>(&X) != NULL ) {
            const ThyraVectorWrapper* x = dynamic_cast<const ThyraVectorWrapper*>(&X);
            AMP_ASSERT(x!=NULL);
            size_t N = d_vecs.size();       // Number of columns of M
            size_t M = x->d_vecs.size();    // Number of columns of X
            AMP_ASSERT(Y->domain()->dim()==(Teuchos::Ordinal)M);
            for (size_t i=0; i<M; i++) {
                Teuchos::RCP<Thyra::VectorBase<double> > Y2 = Y->col(0);
                size_t N2 = static_cast<size_t>( Y2->space()->dim() );
                AMP_ASSERT(d_N_cols==N2);
                for (size_t j=0; j<N; j++) {
                    double dot = d_vecs[j]->dot(x->d_vecs[i]);
                    double y = 0.0;
                    if ( beta!=0.0 )
                        y = beta*Thyra::get_ele(*Y2,d_cols[j]);
                    Thyra::set_ele(d_cols[j],alpha*dot+y,Y2.ptr());
                }
            }        
        } else {
            AMP_ERROR("Not finished");
        }
    } else {
        AMP_ERROR("Unknown case");
    }
}


/****************************************************************
* Protected functions derived from Thyra::MultiVectorBase       *
****************************************************************/
Teuchos::RCP<Thyra::VectorBase<double> > ThyraVectorWrapper::nonconstColImpl(Teuchos::Ordinal j)
{
    std::vector<size_t> cols(1,j);
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vecs(1,d_vecs[j]);
    return Teuchos::RCP<ThyraVectorWrapper>(new ThyraVectorWrapper(vecs,cols,d_N_cols));
}
Teuchos::RCP<const Thyra::MultiVectorBase<double> > ThyraVectorWrapper::contigSubViewImpl(const Teuchos::Range1D &colRng) const
{
    std::vector<size_t> cols;
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vecs;
    for (Teuchos::Ordinal i=colRng.lbound(); i<=colRng.ubound(); i++) {
        cols.push_back( i );
        vecs.push_back( d_vecs[i] );
    }
    return Teuchos::RCP<const ThyraVectorWrapper>(new ThyraVectorWrapper(vecs,cols,d_N_cols));
}
Teuchos::RCP<Thyra::MultiVectorBase<double> > ThyraVectorWrapper::nonconstContigSubViewImpl(const Teuchos::Range1D &colRng)
{
    std::vector<size_t> cols;
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vecs;
    for (Teuchos::Ordinal i=colRng.lbound(); i<=colRng.ubound(); i++) {
        cols.push_back( i );
        vecs.push_back( d_vecs[i] );
    }
    return Teuchos::RCP<ThyraVectorWrapper>(new ThyraVectorWrapper(vecs,cols,d_N_cols));
}
Teuchos::RCP<const Thyra::MultiVectorBase<double> > ThyraVectorWrapper::nonContigSubViewImpl(const Teuchos::ArrayView<const int> &columns) const
{
    std::vector<size_t> cols;
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vecs;
    for (Teuchos::ArrayView<const int>::size_type i=0; i<columns.size(); i++) {
        cols.push_back( columns[i] );
        vecs.push_back( d_vecs[columns[i]] );
    }
    return Teuchos::RCP<const ThyraVectorWrapper>(new ThyraVectorWrapper(vecs,cols,d_N_cols));
}
Teuchos::RCP<Thyra::MultiVectorBase<double> > ThyraVectorWrapper::nonconstNonContigSubViewImpl(const Teuchos::ArrayView<const int> &columns)
{
    std::vector<size_t> cols;
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vecs;
    for (Teuchos::ArrayView<const int>::size_type i=0; i<columns.size(); i++) {
        cols.push_back( columns[i] );
        vecs.push_back( d_vecs[columns[i]] );
    }
    return Teuchos::RCP<ThyraVectorWrapper>(new ThyraVectorWrapper(vecs,cols,d_N_cols));
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
    AMP_ASSERT(d_vecs.size()==1);
    Vector::const_shared_ptr vec = getVec(0);
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
    AMP_ASSERT(d_vecs.size()==1);
    Vector::const_shared_ptr vec = getVec(0);
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
    AMP_ASSERT(d_vecs.size()==1);
    Vector::shared_ptr vec = getVec(0);
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
* applyOpImpl                                                   *
****************************************************************/
void ThyraVectorWrapper::applyOpImpl( const RTOpPack::RTOpT<double> &op, 
        const Teuchos::ArrayView<const Teuchos::Ptr<const Thyra::VectorBase<double> > > &vecs, 
        const Teuchos::ArrayView<const Teuchos::Ptr<Thyra::VectorBase<double> > > &targ_vecs, 
        const Teuchos::Ptr<RTOpPack::ReductTarget> &reduct_obj, 
        const Teuchos::Ordinal global_offset) const
{
    size_t n_blocks = d_vecs[0]->numberOfDataBlocks();
    std::vector<size_t> block_size(n_blocks,0);
    for (size_t i=0; i<n_blocks; i++)
        block_size[i] = d_vecs[0]->sizeOfDataBlock(i);
    // Check that all vectors are compatible and get the pointers
    std::vector<const ThyraVectorWrapper*> vecs_ptr = getConstPtr( block_size, vecs );
    std::vector<ThyraVectorWrapper*> targ_vecs_ptr = getPtr( block_size, targ_vecs );
    // Apply the operation
    Teuchos::RCP<RTOpPack::ReductTarget> reduct_obj2 = op.reduct_obj_create();
    for (size_t k=0; k<d_vecs.size(); k++) {
        for (size_t j=0; j<n_blocks; j++) {
            std::vector<RTOpPack::ConstSubVectorView<double> >  sub_vecs(vecs_ptr.size());
            std::vector<RTOpPack::SubVectorView<double> > targ_sub_vecs(targ_vecs_ptr.size());
            for (size_t i=0; i<vecs_ptr.size(); i++) {
                sub_vecs[i] = RTOpPack::ConstSubVectorView<double>( 
                    Teuchos::ArrayRCP<double>(vecs_ptr[i]->d_vecs[k]->getRawDataBlock<double>(j),0,block_size[j],false) );
            }
            for (size_t i=0; i<targ_vecs_ptr.size(); i++) {
                targ_sub_vecs[i] = RTOpPack::SubVectorView<double>( 
                    Teuchos::ArrayRCP<double>(targ_vecs_ptr[i]->d_vecs[k]->getRawDataBlock<double>(j),0,block_size[j],false) );
            }
            Teuchos::ArrayView<const RTOpPack::ConstSubVectorView<double> > sub_vecs2 = 
                 Teuchos::ArrayView<const RTOpPack::ConstSubVectorView<double> >(sub_vecs);
            Teuchos::ArrayView<const RTOpPack::SubVectorView<double> > targ_sub_vecs2 = 
               Teuchos::ArrayView<const RTOpPack::SubVectorView<double> >(targ_sub_vecs);
            op.apply_op( sub_vecs2, targ_sub_vecs2, reduct_obj2.ptr() );
        }
    }
    // Reduce the result
    if ( reduct_obj.get()!=NULL ) {
        RTOpPack::SPMD_all_reduce<double>( d_comm, op, 1, 
            Teuchos::tuple<const RTOpPack::ReductTarget*>(&*reduct_obj2).getRawPtr(),
            Teuchos::tuple<RTOpPack::ReductTarget*>(&*reduct_obj).getRawPtr() );
    }
    // Change the vector state targ_vecs
    for (size_t i=0; i<targ_vecs_ptr.size(); i++) {
        for (size_t j=0; j<d_vecs.size(); j++) {
            if ( targ_vecs_ptr[i]->d_vecs[j]->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED )
                targ_vecs_ptr[i]->d_vecs[j]->setUpdateStatus( AMP::LinearAlgebra::Vector::LOCAL_CHANGED );
        }
    }
}
void ThyraVectorWrapper::mvMultiReductApplyOpImpl( const RTOpPack::RTOpT<double> &primary_op, 
        const Teuchos::ArrayView< const Teuchos::Ptr< const Thyra::MultiVectorBase<double> > > &multi_vecs, 
        const Teuchos::ArrayView< const Teuchos::Ptr< Thyra::MultiVectorBase<double> > > &targ_multi_vecs, 
        const Teuchos::ArrayView< const Teuchos::Ptr< RTOpPack::ReductTarget > > &reduct_objs, 
        const Teuchos::Ordinal primary_global_offset ) const
{
    size_t n_blocks = d_vecs[0]->numberOfDataBlocks();
    std::vector<size_t> block_size(n_blocks,0);
    for (size_t i=0; i<n_blocks; i++)
        block_size[i] = d_vecs[0]->sizeOfDataBlock(i);
    // Check that all vectors are compatible and get the pointers
    std::vector<const ThyraVectorWrapper*> vecs_ptr = getConstPtr( block_size, multi_vecs );
    std::vector<ThyraVectorWrapper*> targ_vecs_ptr = getPtr( block_size, targ_multi_vecs );
    // Apply the operation
    Teuchos::RCP<RTOpPack::ReductTarget> reduct_obj2 = primary_op.reduct_obj_create();
    for (size_t k=0; k<d_vecs.size(); k++) {
        for (size_t j=0; j<n_blocks; j++) {
            std::vector<RTOpPack::ConstSubVectorView<double> >  sub_vecs(vecs_ptr.size());
            std::vector<RTOpPack::SubVectorView<double> > targ_sub_vecs(targ_vecs_ptr.size());
            for (size_t i=0; i<vecs_ptr.size(); i++) {
                sub_vecs[i] = RTOpPack::ConstSubVectorView<double>( 
                    Teuchos::ArrayRCP<double>(vecs_ptr[i]->d_vecs[k]->getRawDataBlock<double>(j),0,block_size[j],false) );
            }
            for (size_t i=0; i<targ_vecs_ptr.size(); i++) {
                targ_sub_vecs[i] = RTOpPack::SubVectorView<double>( 
                    Teuchos::ArrayRCP<double>(targ_vecs_ptr[i]->d_vecs[k]->getRawDataBlock<double>(j),0,block_size[j],false) );
            }
            Teuchos::ArrayView<const RTOpPack::ConstSubVectorView<double> > sub_vecs2 = 
                 Teuchos::ArrayView<const RTOpPack::ConstSubVectorView<double> >(sub_vecs);
            Teuchos::ArrayView<const RTOpPack::SubVectorView<double> > targ_sub_vecs2 = 
               Teuchos::ArrayView<const RTOpPack::SubVectorView<double> >(targ_sub_vecs);
            primary_op.apply_op( sub_vecs2, targ_sub_vecs2, reduct_obj2.ptr() );
        }
    }
    // Reduce the result
    for (int i=0; i<reduct_objs.size(); i++) {
        RTOpPack::SPMD_all_reduce<double>( d_comm, primary_op, 1, 
            Teuchos::tuple<const RTOpPack::ReductTarget*>(&*reduct_obj2).getRawPtr(),
            Teuchos::tuple<RTOpPack::ReductTarget*>(&*reduct_objs[i]).getRawPtr() );
    }
    // Change the vector state targ_vecs
    for (size_t i=0; i<targ_vecs_ptr.size(); i++) {
        for (size_t j=0; j<d_vecs.size(); j++) {
            if ( targ_vecs_ptr[i]->d_vecs[j]->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED )
                targ_vecs_ptr[i]->d_vecs[j]->setUpdateStatus( AMP::LinearAlgebra::Vector::LOCAL_CHANGED );
        }
    }
}
void ThyraVectorWrapper::mvSingleReductApplyOpImpl( const RTOpPack::RTOpT<double> &primary_op, 
        const RTOpPack::RTOpT<double> &secondary_op, 
        const Teuchos::ArrayView< const Teuchos::Ptr< const Thyra::MultiVectorBase<double> > > &multi_vecs, 
        const Teuchos::ArrayView< const Teuchos::Ptr< Thyra::MultiVectorBase<double> > > &targ_multi_vecs, 
        const Teuchos::Ptr< RTOpPack::ReductTarget > &reduct_obj, 
        const Teuchos::Ordinal primary_global_offset ) const
{
    size_t n_blocks = d_vecs[0]->numberOfDataBlocks();
    std::vector<size_t> block_size(n_blocks,0);
    for (size_t i=0; i<n_blocks; i++)
        block_size[i] = d_vecs[0]->sizeOfDataBlock(i);
    // Check that all vectors are compatible and get the pointers
    std::vector<const ThyraVectorWrapper*> vecs_ptr = getConstPtr( block_size, multi_vecs );
    std::vector<ThyraVectorWrapper*> targ_vecs_ptr = getPtr( block_size, targ_multi_vecs );
    // Apply the operation
    Teuchos::RCP<RTOpPack::ReductTarget> reduct_obj2 = primary_op.reduct_obj_create();
    for (size_t k=0; k<d_vecs.size(); k++) {
        for (size_t j=0; j<n_blocks; j++) {
            std::vector<RTOpPack::ConstSubVectorView<double> >  sub_vecs(vecs_ptr.size());
            std::vector<RTOpPack::SubVectorView<double> > targ_sub_vecs(targ_vecs_ptr.size());
            for (size_t i=0; i<vecs_ptr.size(); i++) {
                sub_vecs[i] = RTOpPack::ConstSubVectorView<double>( 
                    Teuchos::ArrayRCP<double>(vecs_ptr[i]->d_vecs[k]->getRawDataBlock<double>(j),0,block_size[j],false) );
            }
            for (size_t i=0; i<targ_vecs_ptr.size(); i++) {
                targ_sub_vecs[i] = RTOpPack::SubVectorView<double>( 
                    Teuchos::ArrayRCP<double>(targ_vecs_ptr[i]->d_vecs[k]->getRawDataBlock<double>(j),0,block_size[j],false) );
            }
            Teuchos::ArrayView<const RTOpPack::ConstSubVectorView<double> > sub_vecs2 = 
                 Teuchos::ArrayView<const RTOpPack::ConstSubVectorView<double> >(sub_vecs);
            Teuchos::ArrayView<const RTOpPack::SubVectorView<double> > targ_sub_vecs2 = 
               Teuchos::ArrayView<const RTOpPack::SubVectorView<double> >(targ_sub_vecs);
            primary_op.apply_op( sub_vecs2, targ_sub_vecs2, reduct_obj2.ptr() );
        }
    }
    // Reduce the result
    if ( reduct_obj.get()!=NULL ) {
        RTOpPack::SPMD_all_reduce<double>( d_comm, primary_op, 1, 
            Teuchos::tuple<const RTOpPack::ReductTarget*>(&*reduct_obj2).getRawPtr(),
            Teuchos::tuple<RTOpPack::ReductTarget*>(&*reduct_obj).getRawPtr() );
    }
    // Change the vector state targ_vecs
    for (size_t i=0; i<targ_vecs_ptr.size(); i++) {
        for (size_t j=0; j<d_vecs.size(); j++) {
            if ( targ_vecs_ptr[i]->d_vecs[j]->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED )
                targ_vecs_ptr[i]->d_vecs[j]->setUpdateStatus( AMP::LinearAlgebra::Vector::LOCAL_CHANGED );
        }
    }
}


/****************************************************************
* Helper functions                                              *
****************************************************************/
template<class TYPE>
static std::vector<ThyraVectorWrapper*> getPtr( std::vector<size_t> block_size,
    Teuchos::ArrayView<const Teuchos::Ptr<TYPE> > vecs )
{
    int n_vecs = vecs.size();
    std::vector<ThyraVectorWrapper*> ptr(n_vecs,NULL);
    for (int i=0; i<n_vecs; i++) {
        ptr[i] = dynamic_cast<ThyraVectorWrapper*>( vecs[i].getRawPtr() );
        AMP_INSIST(ptr[i]!=NULL,"All vectors used in applyOpImpl must be of the type ThyraVectorWrapper");
        for (size_t j=0; j<ptr[i]->numVecs(); j++) {
            AMP::LinearAlgebra::Vector::const_shared_ptr tmp = ptr[i]->getVec(j);
            AMP_ASSERT(tmp->numberOfDataBlocks()==block_size.size());
            for (size_t k=0; k<block_size.size(); k++)
                AMP_ASSERT(block_size[k]==tmp->sizeOfDataBlock(k));
        }
    }
    return ptr;
}
template<class TYPE>
static std::vector<const ThyraVectorWrapper*> getConstPtr( std::vector<size_t> block_size,
    Teuchos::ArrayView<const Teuchos::Ptr<const TYPE> > vecs )
{
    int n_vecs = vecs.size();
    std::vector<const ThyraVectorWrapper*> ptr(n_vecs,NULL);
    for (int i=0; i<n_vecs; i++) {
        ptr[i] = dynamic_cast<const ThyraVectorWrapper*>( vecs[i].getRawPtr() );
        AMP_INSIST(ptr[i]!=NULL,"All vectors used in applyOpImpl must be of the type ThyraVectorWrapper");
        for (size_t j=0; j<ptr[i]->numVecs(); j++) {
            AMP::LinearAlgebra::Vector::const_shared_ptr tmp = ptr[i]->getVec(j);
            AMP_ASSERT(tmp->numberOfDataBlocks()==block_size.size());
            for (size_t k=0; k<block_size.size(); k++)
                AMP_ASSERT(block_size[k]==tmp->sizeOfDataBlock(k));
        }
    }
    return ptr;
}


}
}

