#include "RNG.h"

namespace AMP {

double RNG::d_SizeTDivisor;
size_t RNG::d_Seed;


void RNG::initialize( size_t seed )
{
    d_Seed         = seed;
    size_t max_int = 0;
    max_int--;
    d_SizeTDivisor = static_cast<double>( max_int );
    double mach_eps, t1, t2;
    t1 = 1;
    t2 = 2;
    while ( t2 > 1. ) {
        mach_eps = t1;
        t1 /= 2.;
        t2 = 1. + t1;
    }
    d_SizeTDivisor *= ( 1. + mach_eps );
}

RNG::RNG( RNGParameters::shared_ptr params ) : d_Params( params )
{
    if ( params->d_WhichSeed == RNGParameters::RNGOptions::USE_GLOBAL_SEED ) {
        srand( static_cast<int>( d_Seed ) + params->d_Rank );
    } else {
        srand( static_cast<int>( params->d_Seed ) + params->d_Rank );
    }
}

RNG::shared_ptr RNG::cloneRNG( size_t new_rank )
{
    AMP_ASSERT( new_rank != d_Params->d_Rank );
    RNGParameters::shared_ptr newParams(
        new RNGParameters( d_Params->d_WhichSeed, new_rank, d_Params->d_Seed ) );
    return RNG::shared_ptr( new RNG( newParams ) );
}
}
