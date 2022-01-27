#include "AMP/vectors/trilinos/epetra/EpetraVectorData.h"
#include "AMP/vectors/data/VectorDataCPU.h"


namespace AMP::LinearAlgebra {

static inline double *getBufferPtr( std::shared_ptr<VectorData> buf )
{
    size_t N_blocks = buf->numberOfDataBlocks();
    if ( N_blocks == 0 )
        return nullptr;
    if ( N_blocks > 1 )
        AMP_ERROR( "More than 1 data block detected" );
    return buf->getRawDataBlock<double>( 0 );
}


EpetraVectorData::EpetraVectorData( std::shared_ptr<EpetraVectorEngineParameters> alias,
                                    Epetra_DataAccess method,
                                    const Epetra_BlockMap &map,
                                    std::shared_ptr<VectorData> bufData,
                                    int localStart,
                                    int localSize,
                                    int globalSize )
    : VectorData( alias->d_CommList ),
      d_epetraVector( method, map, getBufferPtr( bufData ) ),
      d_buf_scope{ bufData }
{
    d_localStart = localStart;
    d_localSize  = localSize;
    d_globalSize = globalSize;
    AMP_ASSERT( d_globalSize < 0x80000000 );
}

std::shared_ptr<EpetraVectorData>
EpetraVectorData::create( std::shared_ptr<EpetraVectorEngineParameters> alias,
                          std::shared_ptr<VectorData> buf )
{
    auto params = std::dynamic_pointer_cast<EpetraVectorEngineParameters>( alias );
    return std::make_shared<EpetraVectorData>( alias,
                                               View,
                                               params->getEpetraMap(),
                                               buf,
                                               params->beginDOF(),
                                               params->getLocalSize(),
                                               params->getGlobalSize() );
}

void *EpetraVectorData::getRawDataBlockAsVoid( size_t i )
{
    if ( i > 1 )
        return nullptr;
    double *p;
    d_epetraVector.ExtractView( &p );
    return p;
}

const void *EpetraVectorData::getRawDataBlockAsVoid( size_t i ) const
{
    if ( i > 1 )
        return nullptr;
    double *p;
    d_epetraVector.ExtractView( &p );
    return p;
}

void EpetraVectorData::setValuesByLocalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ )
        d_epetraVector[indices[i]] = vals[i];
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}

void EpetraVectorData::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    if ( num == 0 )
        return;
    std::vector<int> indices2( num, 0 );
    for ( int i = 0; i < num; i++ )
        indices2[i] = (int) indices[i];
    d_epetraVector.ReplaceGlobalValues( num, const_cast<double *>( vals ), &indices2[0] );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}

void EpetraVectorData::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
    if ( num == 0 )
        return;
    for ( int i = 0; i != num; i++ )
        d_epetraVector[indices[i]] += vals[i];
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}

void EpetraVectorData::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    if ( num == 0 )
        return;
    std::vector<int> indices2( num, 0 );
    for ( int i = 0; i < num; i++ )
        indices2[i] = (int) indices[i];
    d_epetraVector.SumIntoGlobalValues( num, const_cast<double *>( vals ), &indices2[0] );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}

void EpetraVectorData::getValuesByLocalID( int num, size_t *indices, double *vals ) const
{
    if ( num == 0 )
        return;
    double *data;
    d_epetraVector.ExtractView( &data );
    for ( int i = 0; i != num; i++ )
        vals[i] = data[indices[i]];
}

void EpetraVectorData::getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const
{
    if ( num == 0 )
        return;
    double *data;
    d_epetraVector.ExtractView( &data );
    for ( int i = 0; i < num; i++ ) {
        AMP_ASSERT( indices[i] >= d_localStart && indices[i] < d_localStart + d_localSize );
        vals[i] = static_cast<double>( data[indices[i] - d_localStart] );
    }
}

void EpetraVectorData::putRawData( const double *in )
{
    double *p;
    d_epetraVector.ExtractView( &p );
    memcpy( p, in, d_localSize * sizeof( double ) );
}

void EpetraVectorData::copyOutRawData( double *out ) const { d_epetraVector.ExtractCopy( out ); }

std::shared_ptr<VectorData> EpetraVectorData::cloneData() const
{
    auto buffer =
        std::make_shared<VectorDataCPU<double>>( d_localStart, d_localSize, d_globalSize );
    auto params = std::make_shared<EpetraVectorEngineParameters>(
        d_localSize, getComm(), getCommunicationList() );
    auto data = EpetraVectorData::create( params, buffer );
    return data;
}

void EpetraVectorData::swapData( VectorData &rhs )
{
    auto rhs2 = dynamic_cast<EpetraVectorData *>( &rhs );
    AMP_INSIST( rhs2, "Cannot swap with arbitrary VectorData" );
    std::swap( d_epetraVector, rhs2->d_epetraVector );
    std::swap( d_buf_scope, rhs2->d_buf_scope );
    std::swap( d_localStart, rhs2->d_localStart );
    std::swap( d_localSize, rhs2->d_localSize );
    std::swap( d_globalSize, rhs2->d_globalSize );
}


} // namespace AMP::LinearAlgebra
