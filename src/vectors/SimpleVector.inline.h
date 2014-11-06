#include "VectorSelector.h"

namespace AMP {
namespace LinearAlgebra {


inline 
void SimpleVector::copyVector( Vector::const_shared_ptr src_vec )
{
    if ( getLocalSize() != src_vec->getLocalSize() )
        AMP_ERROR( "Mismatched vectors" );
    ConstVectorDataIterator it = src_vec->begin();
    for (size_t i=0; i<getLocalSize(); i++) {
        d_Data[i] = *it;
        ++it;
    }
}


inline
AMP::shared_ptr<ParameterBase> SimpleVector::getParameters ()
{
    AMP_ERROR( "Not implemented" );
    return AMP::shared_ptr<ParameterBase> ();
}


inline
Vector::shared_ptr SimpleVector::cloneVector(const Variable::shared_ptr name) const
{
    return create ( name, d_DOFManager, getCommunicationList() );
}


inline
size_t  SimpleVector::numberOfDataBlocks () const
{
    return 1;
}


inline
size_t  SimpleVector::sizeOfDataBlock ( size_t i ) const
{
    if ( i > 0 )
        return 0;
    return d_Data.size();
}


inline
void SimpleVector::swapVectors(Vector &rhs)
{
    d_Data.swap ( rhs.castTo<SimpleVector>().d_Data );
}


inline
void SimpleVector::aliasVector(Vector &)
{
    AMP_ERROR( "Not implemented" );
}


inline
void SimpleVector::setValuesByLocalID ( int num , size_t *indices , const double *vals )
{
    INCREMENT_COUNT("Virtual");
    for ( int i = 0 ; i != num ; i++ )
        d_Data[indices[i]] = vals[i];
    if ( *d_UpdateState == UNCHANGED )
        *d_UpdateState = LOCAL_CHANGED;
}


inline
void SimpleVector::setLocalValuesByGlobalID ( int num, size_t *indices, const double *vals )
{
    INCREMENT_COUNT("Virtual");
    for ( int i = 0 ; i != num ; i++ ) {
        AMP_ASSERT( indices[i]>=d_startIndex && indices[i]<d_startIndex+d_Data.size() );
        d_Data[indices[i]-d_startIndex] = vals[i];
    }
    if ( *d_UpdateState == UNCHANGED )
        *d_UpdateState = LOCAL_CHANGED;
}


inline
void SimpleVector::addValuesByLocalID ( int num, size_t *indices, const double *vals )
{
    INCREMENT_COUNT("Virtual");
    for ( int i = 0 ; i != num ; i++ )
        d_Data[indices[i]] += vals[i];
    if ( *d_UpdateState == UNCHANGED )
        *d_UpdateState = LOCAL_CHANGED;
}


inline
void SimpleVector::addLocalValuesByGlobalID ( int num, size_t *indices, const double *vals )
{
    INCREMENT_COUNT("Virtual");
    for ( int i = 0 ; i != num ; i++ ) {
        AMP_ASSERT( indices[i]>=d_startIndex && indices[i]<d_startIndex+d_Data.size() );
        d_Data[indices[i]-d_startIndex] += vals[i];
    }
    if ( *d_UpdateState == UNCHANGED )
        *d_UpdateState = LOCAL_CHANGED;
}


inline
void SimpleVector::getLocalValuesByGlobalID ( int num, size_t *indices, double *vals ) const
{
    INCREMENT_COUNT("Virtual");
    for ( int i = 0 ; i != num ; i++ ) {
        AMP_ASSERT( indices[i]>=d_startIndex && indices[i]<d_startIndex+d_Data.size() );
        vals[i] = d_Data[indices[i]-d_startIndex];
    }
}


inline
void SimpleVector::assemble()
{
    AMP_ERROR( "Not implemented" );
}


inline
size_t SimpleVector::getLocalSize() const
{
    return d_Data.size();
}


inline
size_t SimpleVector::getGlobalSize() const
{
    return d_globalSize;
}


inline
void *SimpleVector::getRawDataBlockAsVoid ( size_t i )
{
    if ( i != 0 )
    {
      return 0;
    }
    return &(d_Data[0]);
}


inline
const void *SimpleVector::getRawDataBlockAsVoid ( size_t i ) const
{
    if ( i != 0 )
    {
      return 0;
    }
    return &(d_Data[0]);
}


inline
double &SimpleVector::operator[] ( size_t i )
{
    AMP_ASSERT ( i < d_Data.size() );
    return d_Data[i];
}


inline
double  SimpleVector::operator[] ( size_t i ) const 
{
    AMP_ASSERT ( i < d_Data.size() );
    return d_Data[i];
}


inline
void SimpleVector::resize ( size_t i )
{
    d_Data.resize ( i );
}


}
}
