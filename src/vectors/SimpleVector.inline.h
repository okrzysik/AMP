#include "VectorSelector.h"

namespace AMP {
namespace LinearAlgebra {

template <typename T>
inline Vector::shared_ptr SimpleVector<T>::cloneVector( const Variable::shared_ptr name ) const
{
    return create( name, d_DOFManager, getCommunicationList() );
}

template <typename T>
inline size_t SimpleVector<T>::numberOfDataBlocks() const
{
    return 1;
}

template <typename T>
inline size_t SimpleVector<T>::sizeOfDataBlock( size_t i ) const
{
    if ( i > 0 ) return 0;
    return d_Data.size();
}

template <typename T>
inline void SimpleVector<T>::swapVectors( Vector &rhs )
{
    d_Data.swap( rhs.castTo<SimpleVector>().d_Data );
}

template <typename T>
inline void SimpleVector<T>::aliasVector( Vector & )
{
    AMP_ERROR( "Not implemented" );
}

template <typename T>
inline void SimpleVector<T>::setValuesByLocalID( int num, size_t *indices, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
    for ( int i = 0; i != num; i++ ) d_Data[indices[i]] = static_cast<T>( vals[i] );
    if ( *d_UpdateState == UNCHANGED ) *d_UpdateState   = LOCAL_CHANGED;
}

template <typename T>
inline void
SimpleVector<T>::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
    for ( int i = 0; i != num; i++ ) {
        AMP_ASSERT( indices[i] >= d_startIndex && indices[i] < d_startIndex + d_Data.size() );
        d_Data[indices[i] - d_startIndex] = static_cast<T>( vals[i] );
    }
    if ( *d_UpdateState == UNCHANGED ) *d_UpdateState = LOCAL_CHANGED;
}

template <typename T>
inline void SimpleVector<T>::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
    for ( int i = 0; i != num; i++ ) d_Data[indices[i]] += static_cast<T>( vals[i] );
    if ( *d_UpdateState == UNCHANGED ) *d_UpdateState = LOCAL_CHANGED;
}

template <typename T>
inline void
SimpleVector<T>::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
    for ( int i = 0; i != num; i++ ) {
        AMP_ASSERT( indices[i] >= d_startIndex && indices[i] < d_startIndex + d_Data.size() );
        d_Data[indices[i] - d_startIndex] += static_cast<T>( vals[i] );
    }
    if ( *d_UpdateState == UNCHANGED ) *d_UpdateState = LOCAL_CHANGED;
}

template <typename T>
inline void
SimpleVector<T>::getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const
{
    INCREMENT_COUNT( "Virtual" );
    for ( int i = 0; i != num; i++ ) {
        AMP_ASSERT( indices[i] >= d_startIndex && indices[i] < d_startIndex + d_Data.size() );
        vals[i] = static_cast<double>( d_Data[indices[i] - d_startIndex] );
    }
}

template <typename T>
inline void SimpleVector<T>::assemble()
{
    AMP_ERROR( "Not implemented" );
}

template <typename T>
inline size_t SimpleVector<T>::getLocalSize() const
{
    return d_Data.size();
}

template <typename T>
inline size_t SimpleVector<T>::getGlobalSize() const
{
    return d_globalSize;
}

template <typename T>
inline void *SimpleVector<T>::getRawDataBlockAsVoid( size_t i )
{
    if ( i != 0 ) {
        return 0;
    }
    return d_Data.data();
}

template <typename T>
inline const void *SimpleVector<T>::getRawDataBlockAsVoid( size_t i ) const
{
    if ( i != 0 ) {
        return 0;
    }
    return d_Data.data();
}

template <typename T>
inline T &SimpleVector<T>::operator[]( size_t i )
{
    AMP_ASSERT( i < d_Data.size() );
    return d_Data[i];
}

template <typename T>
inline T SimpleVector<T>::operator[]( size_t i ) const
{
    AMP_ASSERT( i < d_Data.size() );
    return d_Data[i];
}

template <typename T>
inline void SimpleVector<T>::resize( size_t i )
{
    d_Data.resize( i );
}
}
}
