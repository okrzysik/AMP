#include "AMP/vectors/trilinos/epetra/EpetraVectorData.h"
#include "AMP/vectors/data/VectorDataCPU.h"


namespace AMP {
namespace LinearAlgebra {

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
  : VectorData(alias),
    d_epetraVector( method, map, getBufferPtr(bufData) ),
    d_buf_scope{ bufData },
    d_iLocalStart{ localStart },
    d_iLocalSize{ localSize },
    d_iGlobalSize{ globalSize }
{
  
}

std::shared_ptr<EpetraVectorData> EpetraVectorData::create( std::shared_ptr<EpetraVectorEngineParameters> alias,
							    std::shared_ptr<VectorData> buf)
{
  return std::make_shared<EpetraVectorData>(alias,
					    View,
					    std::dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->getEpetraMap(),
					    buf,
					    std::dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->beginDOF(),
					    std::dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->getLocalSize(),
					    std::dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->getGlobalSize());
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
}

void EpetraVectorData::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    if ( num == 0 )
        return;
    AMP_ASSERT( getGlobalSize() < 0x80000000 );
    std::vector<int> indices2( num, 0 );
    for ( int i = 0; i < num; i++ )
        indices2[i] = (int) indices[i];
    d_epetraVector.ReplaceGlobalValues( num, const_cast<double *>( vals ), &indices2[0] );
}

void EpetraVectorData::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
    if ( num == 0 )
        return;
    for ( int i = 0; i != num; i++ )
        d_epetraVector[indices[i]] += vals[i];
}

void EpetraVectorData::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    if ( num == 0 )
        return;
    AMP_ASSERT( getGlobalSize() < 0x80000000 );
    std::vector<int> indices2( num, 0 );
    for ( int i = 0; i < num; i++ )
        indices2[i] = (int) indices[i];
    d_epetraVector.SumIntoGlobalValues( num, const_cast<double *>( vals ), &indices2[0] );
}

void EpetraVectorData::getValuesByLocalID( int, size_t *, double * ) const
{
    AMP_ERROR( "EpetraVectorData::getValuesByLocalID() This shouldn't be called" );
}

void EpetraVectorData::getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const
{
    //    AMP_ERROR( "EpetraVectorData::getLocalValuesByGlobalID() This shouldn't be called" );
    if ( num == 0 )
        return;
    double *data;
    d_epetraVector.ExtractView( &data );
    for ( int i = 0; i < num; i++ ) {
        AMP_ASSERT( (int64_t) indices[i] >= d_iLocalStart &&
                    (int64_t) indices[i] < d_iLocalStart + d_iLocalSize );
        vals[i] = static_cast<double>( data[indices[i] - d_iLocalStart] );
    }
}

void EpetraVectorData::putRawData( const double *in )
{
    double *p;
    d_epetraVector.ExtractView( &p );
    size_t N = getLocalSize();
    memcpy( p, in, N * sizeof( double ) );
}

void EpetraVectorData::copyOutRawData( double *out ) const { d_epetraVector.ExtractCopy( out ); }

std::shared_ptr<VectorData> EpetraVectorData::cloneData() const
{
    auto buffer = std::make_shared<VectorDataCPU<double>>(
        getLocalStartID(), getLocalSize(), getGlobalSize() );
    return buffer;
}


} // namespace LinearAlgebra
} // namespace AMP
