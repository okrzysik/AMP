#ifndef included_FunctionTable_hpp
#define included_FunctionTable_hpp

#include "utils/FunctionTable.h"
#include "utils/Utilities.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <random>


namespace AMP {


/********************************************************
*  Random number initialization                         *
********************************************************/
template <class TYPE, class FUN>
void FunctionTable::rand( Array<TYPE, FUN> &x )
{
    FunctionTable::rand<TYPE>( x.length(), x.data() );
}
template <>
inline void FunctionTable::rand<double>( size_t N, double *x )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<> dis( 0, 1 );
    for ( size_t i = 0; i < N; i++ )
        x[i]       = dis( gen );
}
template <>
inline void FunctionTable::rand<float>( size_t N, float *x )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<> dis( 0, 1 );
    for ( size_t i = 0; i < N; i++ )
        x[i]       = dis( gen );
}
template <>
inline void FunctionTable::rand<int>( size_t N, int *x )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_int_distribution<> dis;
    for ( size_t i = 0; i < N; i++ )
        x[i]       = dis( gen );
}


/********************************************************
*  Reduction                                            *
********************************************************/
template <class TYPE, class FUN, typename LAMBDA>
inline TYPE FunctionTable::reduce( LAMBDA &op, const Array<TYPE, FUN> &A )
{
    if ( A.length() == 0 )
        return TYPE();
    const TYPE *x = A.data();
    TYPE y        = x[0];
    for ( size_t i = 1; i < A.length(); i++ )
        y = op( x[i], y );
    return y;
}


/********************************************************
*  Unary transformation                                 *
********************************************************/
template <class TYPE, class FUN, typename LAMBDA>
inline void FunctionTable::transform( LAMBDA &fun, const Array<TYPE, FUN> &x, Array<TYPE, FUN> &y )
{
    y.resize( x.size() );
    for ( size_t i = 0; i < x.length(); i++ )
        y( i )     = fun( x( i ) );
}
template <class TYPE, class FUN, typename LAMBDA>
inline void FunctionTable::transform( LAMBDA &fun,
                                      const Array<TYPE, FUN> &x,
                                      const Array<TYPE, FUN> &y,
                                      Array<TYPE, FUN> &z )
{
    ARRAY_INSIST( x.sizeMatch( y ), "Sizes of x and y do not match" );
    z.resize( x.size() );
    for ( size_t i = 0; i < x.length(); i++ )
        z( i )     = fun( x( i ), y( i ) );
}


/********************************************************
*  Multiply two arrays                                  *
********************************************************/
template <class TYPE, class FUN>
void FunctionTable::multiply( const Array<TYPE, FUN> &a,
                              const Array<TYPE, FUN> &b,
                              Array<TYPE, FUN> &c )
{
    if ( a.d_ndim <= 2 && b.d_ndim <= 2 ) {
        ARRAY_INSIST( a.size( 1 ) == b.size( 0 ), "Inner dimensions must match" );
        c.resize( a.size( 0 ), b.size( 1 ) );
        c.fill( 0 );
        for ( size_t k = 0; k < b.size( 1 ); k++ ) {
            for ( size_t j = 0; j < a.size( 1 ); j++ ) {
                for ( size_t i = 0; i < a.size( 0 ); i++ ) {
                    c( i, k ) += a( i, j ) * b( j, k );
                }
            }
        }
    } else {
        ARRAY_ERROR( "Not finished yet" );
    }
}


} // namespace AMP

#endif
