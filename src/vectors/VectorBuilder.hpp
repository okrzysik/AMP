#ifndef included_AMP_VectorBuider_hpp
#define included_AMP_VectorBuider_hpp

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/data/ArrayVectorData.h"
#ifdef USE_DEVICE
    #include "AMP/utils/device/GPUFunctionTable.h"
    #include "AMP/vectors/operations/VectorOperationsDevice.h"
#endif
#include "AMP/utils/memory.h"

#include "math.h"


namespace AMP::LinearAlgebra {


/****************************************************************
 * SimpleVector                                                  *
 ****************************************************************/
template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr createSimpleVector( size_t localSize, const std::string &name )
{
    AMP_MPI comm( AMP_COMM_SELF );
    auto var = std::make_shared<Variable>( name );
    return createSimpleVector<TYPE, OPS, DATA>( localSize, var, comm );
}

template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr createSimpleVector( size_t localSize, std::shared_ptr<Variable> var )
{
    AMP_MPI comm( AMP_COMM_SELF );
    return createSimpleVector<TYPE, OPS, DATA>( localSize, var, comm );
}
template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr
createSimpleVector( size_t localSize, std::shared_ptr<Variable> var, AMP_MPI comm )
{
    auto DOFs = std::make_shared<AMP::Discretization::DOFManager>( localSize, comm );
    auto ops  = std::make_shared<OPS>();
    auto data =
        std::make_shared<DATA>( DOFs->beginDOF(), DOFs->numLocalDOF(), DOFs->numGlobalDOF() );
    data->setCommunicationList( std::make_shared<CommunicationList>( localSize, comm ) );
    return std::make_shared<Vector>( data, ops, var, DOFs );
}
template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr createSimpleVector( std::shared_ptr<Variable> var,
                                       std::shared_ptr<AMP::Discretization::DOFManager> DOFs,
                                       std::shared_ptr<CommunicationList> commlist )
{
    auto ops = std::make_shared<OPS>();
    auto data =
        std::make_shared<DATA>( DOFs->beginDOF(), DOFs->numLocalDOF(), DOFs->numGlobalDOF() );
    data->setCommunicationList( commlist );
    return std::make_shared<Vector>( data, ops, var, DOFs );
}


/****************************************************************
 * ArrayVector                                                  *
 ****************************************************************/
template<typename T, typename FUN, typename Allocator>
Vector::shared_ptr createArrayVector( const ArraySize &localSize, const std::string &name )
{
    auto var = std::make_shared<Variable>( name );
    return createArrayVector<T, FUN, Allocator>( localSize, { 0, 0, 0 }, AMP_COMM_SELF, var );
}
template<typename T, typename FUN, typename Allocator>
Vector::shared_ptr createArrayVector( const ArraySize &localSize, std::shared_ptr<Variable> var )
{
    return createArrayVector<T, FUN, Allocator>( localSize, { 0, 0, 0 }, AMP_COMM_SELF, var );
}
template<typename T, typename FUN, typename Allocator>
Vector::shared_ptr createArrayVector( const ArraySize &localSize,
                                      const ArraySize &blockIndex,
                                      const AMP_MPI &comm,
                                      std::shared_ptr<Variable> var )
{
    size_t N  = localSize.length();
    auto ops  = std::make_shared<VectorOperationsDefault<T>>();
    auto data = ArrayVectorData<T, FUN, Allocator>::create( localSize, blockIndex, comm );
    auto DOFs = std::make_shared<AMP::Discretization::DOFManager>( N, comm );
    return std::make_shared<Vector>( data, ops, var, DOFs );
}


/****************************************************************
 * Vector view based on ArrayVector                             *
 ****************************************************************/
template<typename T>
Vector::shared_ptr createVectorAdaptor( const std::string &name,
                                        std::shared_ptr<AMP::Discretization::DOFManager> DOFs,
                                        T *data )
{
    auto var              = std::make_shared<Variable>( name );
    auto params           = std::make_shared<CommunicationListParameters>();
    params->d_comm        = DOFs->getComm();
    params->d_localsize   = DOFs->numLocalDOF();
    params->d_remote_DOFs = DOFs->getRemoteDOFs();
    auto commList         = std::make_shared<CommunicationList>( params );

    auto memType = AMP::Utilities::getMemoryType( data );

    std::shared_ptr<VectorData> vecData;
    std::shared_ptr<VectorOperations> vecOps;
    if ( memType <= AMP::Utilities::MemoryType::host ) {
        vecOps  = std::make_shared<VectorOperationsDefault<T>>();
        vecData = ArrayVectorData<T>::create( DOFs->numLocalDOF(), commList, data );
    } else if ( memType == AMP::Utilities::MemoryType::managed ) {
#ifdef USE_DEVICE
        vecOps = std::make_shared<VectorOperationsDevice<T>>();
#endif
        vecData = ArrayVectorData<T, AMP::GPUFunctionTable, AMP::ManagedAllocator<T>>::create(
            DOFs->numLocalDOF(), commList, data );
    } else if ( memType == AMP::Utilities::MemoryType::device ) {
#ifdef USE_DEVICE
        vecOps = std::make_shared<VectorOperationsDevice<T>>();
#endif
        vecData = ArrayVectorData<T, AMP::GPUFunctionTable, AMP::DeviceAllocator<T>>::create(
            DOFs->numLocalDOF(), commList, data );
    } else {
        AMP_ERROR( "Unknown memory location specified for data" );
    }

    auto vec = std::make_shared<Vector>( vecData, vecOps, var, DOFs );
    vec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    return vec;
}


} // namespace AMP::LinearAlgebra

#endif
