#include "math.h"

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"
#include "ArrayVector.h"

namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
template<typename T, typename FUN, typename Allocator>
ArrayVector<T, FUN, Allocator>::ArrayVector() : Vector()
{
    d_VectorOps = std::make_shared<VectorOperationsDefault<T>>();
}

template<typename T, typename FUN, typename Allocator>
ArrayVector<T, FUN, Allocator>::ArrayVector( std::shared_ptr<ArrayVectorData<T, FUN, Allocator>> data ) : Vector(), d_VectorDataSP{ data }
{
    d_VectorData = d_VectorDataSP.get();
    d_VectorOps = std::make_shared<VectorOperationsDefault<T>>();
}

template<typename T, typename FUN, typename Allocator>
Vector::shared_ptr ArrayVector<T, FUN, Allocator>::create( const std::vector<size_t> &localSize,
                                                           Variable::shared_ptr var )
{
    auto data = ArrayVectorData<T, FUN, Allocator>::create(localSize);
    auto vdata = std::dynamic_pointer_cast<ArrayVectorData<T, FUN, Allocator>>(data);
    auto retVal = std::make_shared<ArrayVector<T, FUN, Allocator>>(vdata);
    retVal->setVariable( var );
    const auto N         = vdata->getArray().length();
    AMP_MPI comm( AMP_COMM_SELF );
    auto DOFs            = std::make_shared<AMP::Discretization::DOFManager>( N, data->getComm() );
    retVal->d_DOFManager = DOFs;
    return retVal;
}

template<typename T, typename FUN, typename Allocator>
Vector::shared_ptr ArrayVector<T, FUN, Allocator>::create( const std::vector<size_t> &localSize,
                                                           Variable::shared_ptr var,
                                                           AMP_MPI comm )
{
    auto data = ArrayVectorData<T, FUN, Allocator>::create(localSize,  comm);
    auto vdata = std::dynamic_pointer_cast<ArrayVectorData<T, FUN, Allocator>>(data);
    auto retVal = std::make_shared<ArrayVector<T, FUN, Allocator>>(vdata);
    retVal->setVariable( var );
    const auto N         = vdata->getArray().length();
    auto DOFs            = std::make_shared<AMP::Discretization::DOFManager>( N, comm );
    retVal->d_DOFManager = DOFs;
    return retVal;
}

template<typename T, typename FUN, typename Allocator>
Vector::shared_ptr
ArrayVector<T, FUN, Allocator>::create( Variable::shared_ptr var,
                                        AMP::Discretization::DOFManager::shared_ptr DOFs,
                                        AMP::LinearAlgebra::CommunicationList::shared_ptr commlist )
{
    auto data = ArrayVectorData<T, FUN, Allocator>::create(commlist);
    auto vdata = std::dynamic_pointer_cast<ArrayVectorData<T, FUN, Allocator>>(data);
    auto retVal = std::make_shared<ArrayVector<T, FUN, Allocator>>(vdata);
    retVal->setVariable( var );
    retVal->d_DOFManager = DOFs;
    AMP_ERROR( "This routine is not complete" );
    return retVal;
}

template<typename T, typename FUN, typename Allocator>
inline Vector::shared_ptr
ArrayVector<T, FUN, Allocator>::cloneVector( const Variable::shared_ptr name ) const
{
    auto vdata = std::dynamic_pointer_cast<ArrayVectorData<T, FUN, Allocator>>(d_VectorDataSP);
    const auto &array = vdata->getArray();
    std::vector<size_t> size( array.size().begin(), array.size().end() );
    return create( size, name, this->getComm() );
}

template<typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::swapVectors( Vector &rhs )
{
    // get internal arrays
    auto vdata = std::dynamic_pointer_cast<ArrayVectorData<T, FUN, Allocator>>(d_VectorDataSP);
    auto &internalArray = vdata->getArray();

    auto &otherArray    = dynamic_cast<ArrayVectorData<T, FUN, Allocator>*>( rhs.getVectorData() )->getArray();
    // reset views
    internalArray.swap( otherArray );
}

template<typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::aliasVector( Vector & )
{
    AMP_ERROR( "Not implemented" );
}

} // namespace LinearAlgebra
} // namespace AMP
