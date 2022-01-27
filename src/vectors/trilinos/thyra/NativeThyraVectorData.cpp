#include "AMP/vectors/trilinos/thyra/NativeThyraVectorData.h"
#include "AMP/vectors/trilinos/thyra/NativeThyraVectorOperations.h"
#include "AMP/vectors/trilinos/thyra/ThyraVectorWrapper.h"


// Trilinos includes
DISABLE_WARNINGS
#include "Thyra_DefaultSpmdVector_def.hpp"
#include "Thyra_VectorStdOps_def.hpp"
#include "Trilinos_version.h"
ENABLE_WARNINGS


namespace AMP::LinearAlgebra {


/************************************************************************
 * Constructors                                                          *
 ************************************************************************/
NativeThyraVectorData::NativeThyraVectorData( Teuchos::RCP<Thyra::VectorBase<double>> vec,
                                              size_t localsize,
                                              AMP_MPI comm )
    : VectorData()
{
    size_t dim = vec->space()->dim();
    AMP_ASSERT( comm.sumReduce( localsize ) == dim );
    auto communicationListParams         = std::make_shared<CommunicationListParameters>();
    communicationListParams->d_comm      = comm;
    communicationListParams->d_localsize = localsize;
    d_CommList   = std::make_shared<CommunicationList>( communicationListParams );
    d_thyraVec   = vec;
    d_localSize  = localsize;
    d_globalSize = d_thyraVec->space()->dim();
    d_localStart = d_CommList->getStartGID();
}


/************************************************************************
 * Destructor                                                            *
 ************************************************************************/
NativeThyraVectorData::~NativeThyraVectorData() = default;

/************************************************************************
 * Vector functions                                                      *
 ************************************************************************/
std::shared_ptr<VectorData> NativeThyraVectorData::cloneData() const
{
    return std::make_shared<NativeThyraVectorData>( d_thyraVec->clone_v(), d_localSize, getComm() );
}

void NativeThyraVectorData::putRawData( const double *in )
{
    size_t i = 0;
    for ( size_t b = 0; b < numberOfDataBlocks(); b++ ) {
        auto *data = reinterpret_cast<double *>( getRawDataBlockAsVoid( b ) );
        for ( size_t j = 0; j < sizeOfDataBlock( b ); j++, i++ )
            data[j] = in[i];
    }
}


void NativeThyraVectorData::copyOutRawData( double *out ) const
{
    size_t i = 0;
    for ( size_t b = 0; b < numberOfDataBlocks(); b++ ) {
        const auto *data = reinterpret_cast<const double *>( getRawDataBlockAsVoid( b ) );
        for ( size_t j = 0; j < sizeOfDataBlock( b ); j++, i++ )
            out[i] = data[j];
    }
}


void *NativeThyraVectorData::getRawDataBlockAsVoid( size_t i )
{
    Thyra::VectorBase<double> *ptr = d_thyraVec.get();
    auto *spmdVector               = dynamic_cast<Thyra::DefaultSpmdVector<double> *>( ptr );
    if ( spmdVector != nullptr ) {
        if ( i != 0 )
            AMP_ERROR( "Invalid block" );
        return spmdVector->getPtr();
    }
    auto *wrapperVector = dynamic_cast<ThyraVectorWrapper *>( ptr );
    if ( wrapperVector != nullptr ) {
        AMP_INSIST( wrapperVector->numVecs() == 1,
                    "Not ready for dealing with multiple copies of the vector yet" );
        return wrapperVector->getVec( 0 )->getRawDataBlock<double>( i );
    }
    AMP_ERROR( "not finished" );
    return nullptr;
}


const void *NativeThyraVectorData::getRawDataBlockAsVoid( size_t i ) const
{
    const Thyra::VectorBase<double> *ptr = d_thyraVec.get();
    const auto *spmdVector = dynamic_cast<const Thyra::DefaultSpmdVector<double> *>( ptr );
    if ( spmdVector != nullptr ) {
        if ( i != 0 )
            AMP_ERROR( "Invalid block" );
        return spmdVector->getPtr();
    }
    const auto *wrapperVector = dynamic_cast<const ThyraVectorWrapper *>( ptr );
    if ( wrapperVector != nullptr ) {
        AMP_INSIST( wrapperVector->numVecs() == 1,
                    "Not ready for dealing with multiple copies of the vector yet" );
        return wrapperVector->getVec( 0 )->getRawDataBlock<double>( i );
    }
    return nullptr;
}


size_t NativeThyraVectorData::numberOfDataBlocks() const
{
    const Thyra::VectorBase<double> *ptr = d_thyraVec.get();
    if ( dynamic_cast<const Thyra::DefaultSpmdVector<double> *>( ptr ) != nullptr )
        return 1;
    const auto *wrapperVector = dynamic_cast<const ThyraVectorWrapper *>( ptr );
    if ( wrapperVector != nullptr )
        return wrapperVector->getVec( 0 )->numberOfDataBlocks();
    AMP_ERROR( "not finished" );
    return 1;
}


size_t NativeThyraVectorData::sizeOfDataBlock( size_t i ) const
{
    const Thyra::VectorBase<double> *ptr = d_thyraVec.get();
    const auto *spmdVector = dynamic_cast<const Thyra::DefaultSpmdVector<double> *>( ptr );
    if ( spmdVector != nullptr )
        return d_localSize;
    const auto *wrapperVector = dynamic_cast<const ThyraVectorWrapper *>( ptr );
    if ( wrapperVector != nullptr ) {
        AMP_INSIST( wrapperVector->numVecs() == 1,
                    "Not ready for dealing with multiple copies of the vector yet" );
        return wrapperVector->getVec( 0 )->sizeOfDataBlock( i );
    }
    AMP_ERROR( "not finished" );
    return d_localSize;
}

Teuchos::RCP<const Thyra::VectorBase<double>>
NativeThyraVectorData::getThyraVec( std::shared_ptr<const VectorData> vec )
{
    auto vec2 = std::dynamic_pointer_cast<const ThyraVector>( vec );
    AMP_ASSERT( vec2 != nullptr );
    return vec2->getVec();
}

Teuchos::RCP<const Thyra::VectorBase<double>>
NativeThyraVectorData::getThyraVec( const VectorData &v )
{
    auto vec2 = dynamic_cast<const ThyraVector *>( &v );
    AMP_ASSERT( vec2 != nullptr );
    return vec2->getVec();
}

Teuchos::RCP<Thyra::VectorBase<double>> NativeThyraVectorData::getThyraVec( VectorData &v )
{
    auto vec2 = dynamic_cast<ThyraVector *>( &v );
    AMP_ASSERT( vec2 != nullptr );
    return vec2->getVec();
}

void NativeThyraVectorData::setValuesByLocalID( int num, size_t *indices, const double *vals )
{
    NULL_USE( num );
    NULL_USE( indices );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


void NativeThyraVectorData::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    NULL_USE( num );
    NULL_USE( indices );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


void NativeThyraVectorData::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
    NULL_USE( num );
    NULL_USE( indices );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


void NativeThyraVectorData::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    NULL_USE( num );
    NULL_USE( indices );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


void NativeThyraVectorData::getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    NULL_USE( numVals );
    NULL_USE( ndx );
    NULL_USE( vals );
    AMP_ERROR( "not implemented" );
}


void NativeThyraVectorData::swapData( VectorData &rhs )
{
    auto rhs2 = dynamic_cast<NativeThyraVectorData *>( &rhs );
    AMP_INSIST( rhs2, "Cannot swap with arbitrary VectorData" );
    std::swap( d_localSize, rhs2->d_localSize );
    std::swap( d_globalSize, rhs2->d_globalSize );
    std::swap( d_localStart, rhs2->d_localStart );
    std::swap( d_thyraVec, rhs2->d_thyraVec );
}


} // namespace AMP::LinearAlgebra
