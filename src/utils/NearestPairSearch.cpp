#include "AMP/utils/NearestPairSearch.h"
#include "AMP/utils/MeshPoint.h"
#include "AMP/utils/NearestPairSearch.hpp"

namespace AMP {


// Calculate the closest pair of points in a list
std::pair<int, int> find_min_dist( const std::vector<AMP::Mesh::MeshPoint<double>> &x )
{
    std::pair<int, int> index( 0, 0 );
    if ( x.empty() )
        return index;
    int ndim = x[0].ndim();
    auto x2  = new double[ndim * x.size()];
    for ( size_t i = 0; i < x.size(); i++ ) {
        for ( int d = 0; d < ndim; d++ )
            x2[d + i * ndim] = x[i][d];
    }
    if ( ndim == 1 )
        index = find_min_dist<1, double>( x.size(), x2 );
    else if ( ndim == 2 )
        index = find_min_dist<2, double>( x.size(), x2 );
    else if ( ndim == 3 )
        index = find_min_dist<3, double>( x.size(), x2 );
    else
        AMP_ERROR( "Not programmed" );
    return index;
}


} // namespace AMP
