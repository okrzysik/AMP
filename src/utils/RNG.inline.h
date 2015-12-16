
#include "utils/Utilities.h"

namespace AMP {

inline RNGParameters::RNGParameters( RNGOptions o, size_t rank, size_t seed )
{
    d_Seed      = seed;
    d_Rank      = rank;
    d_WhichSeed = o;
}

inline void RNG::fillBuffer( void *buf, size_t len )
{
    unsigned char *alias = static_cast<unsigned char *>( buf );
    for ( size_t i = 0; i != len; i++ ) {
        alias[i] = static_cast<unsigned char>( rand() % 256 );
    }
}

inline int RNG::nextInt( int low, int high )
{
    AMP_ASSERT( high > low );
    int buffer;
    fillBuffer( static_cast<void *>( &buffer ), sizeof( int ) );
    return ( abs( buffer ) % ( high - low ) ) + low;
}

inline double RNG::nextDouble( double low, double high )
{
    AMP_ASSERT( high > low );
    return low + ( high - low ) * static_cast<double>( nextInt( 0, 1073741824 ) ) / 1073741824.0;
}

inline RandomVariable<double>::RandomVariable( type low, type high, RNG::shared_ptr r )
    : d_Low( low ), d_High( high ), d_RNG( r )
{
    AMP_ASSERT( high > low );
}


inline RandomVariable<double>::operator type() { return d_RNG->nextDouble( d_Low, d_High ); }

inline RandomVariable<float>::RandomVariable( type low, type high, RNG::shared_ptr r )
    : d_Low( low ), d_High( high ), d_RNG( r )
{
    AMP_ASSERT( high > low );
}


inline RandomVariable<float>::operator type()
{
    return static_cast<float>( d_RNG->nextDouble( d_Low, d_High ) );
}
}
