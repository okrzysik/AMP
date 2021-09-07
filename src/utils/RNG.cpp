#include "AMP/utils/RNG.h"


namespace AMP {

double RNG::d_SizeTDivisor;
size_t RNG::d_Seed;


void RNG::initialize( size_t seed )
{
    d_Seed         = seed;
    size_t max_int = 0;
    max_int--;
    d_SizeTDivisor  = static_cast<double>( max_int );
    double mach_eps = 0;
    double t1       = 1;
    double t2       = 2;
    while ( t2 > 1. ) {
        mach_eps = t1;
        t1 /= 2.;
        t2 = 1. + t1;
    }
    d_SizeTDivisor *= ( 1. + mach_eps );
}

RNG::RNG( std::shared_ptr<RNGParameters> params ) : d_Params( params )
{
    if ( params->d_WhichSeed == RNGParameters::RNGOptions::USE_GLOBAL_SEED ) {
        srand( static_cast<int>( d_Seed ) + params->d_Rank );
    } else {
        srand( static_cast<int>( params->d_Seed ) + params->d_Rank );
    }
}

std::shared_ptr<RNG> RNG::cloneRNG( size_t new_rank )
{
    AMP_ASSERT( new_rank != d_Params->d_Rank );
    auto newParams =
        std::make_shared<RNGParameters>( d_Params->d_WhichSeed, new_rank, d_Params->d_Seed );
    return std::make_shared<RNG>( newParams );
}
} // namespace AMP
