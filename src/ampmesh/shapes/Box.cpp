#include "ampmesh/shapes/Box.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Geometry {


/********************************************************
* Compute the distance to the object                    *
********************************************************/
double Box::distance( const std::initializer_list<double>& position,
    const std::initializer_list<double>& direction ) const
{
    size_t ndim = d_range.size()/2;
    AMP_ASSERT( ndim==position.size() && ndim==direction.size() && ndim<8 );
    double pos[8], dir[8];
    size_t d = 0;
    for ( auto tmp : position )
        pos[d++] = tmp;
    d = 0;
    for ( auto tmp : direction )
        dir[d++] = tmp;
    double dist = std::numeric_limits<double>::quiet_NaN();
    if ( Box::inside(position) ) {
        dist = std::numeric_limits<double>::infinity();
        for (size_t d=0; d<ndim; d++) {
            double d1 = (d_range[2*d+0]-pos[d])/dir[d];
            double d2 = (d_range[2*d+1]-pos[d])/dir[d];
            dist = d1>=0 ? std::min(dist,d1):dist;
            dist = d2>=0 ? std::min(dist,d2):dist;
        }
        dist = -dist;
    } else {
        for (size_t d=0; d<ndim; d++) {
            double d1 = (d_range[2*d+0]-pos[d])/dir[d];
            double d2 = (d_range[2*d+1]-pos[d])/dir[d];
            double dist2 = std::numeric_limits<double>::infinity();
            dist2 = d1>=0 ? std::min(dist2,d1):dist2;
            dist2 = d2>=0 ? std::min(dist2,d2):dist2;
            dist = std::max(dist,dist2);
        }
    }
    return dist;
}


/********************************************************
* Check if the ray is inside the geometry               *
********************************************************/
bool Box::inside( const std::initializer_list<double>& position ) const
{
    size_t ndim = d_range.size()/2;
    AMP_ASSERT( ndim==position.size() && ndim<8 );
    double pos[8];
    size_t d = 0;
    for ( auto tmp : position )
        pos[d++] = tmp;
    bool inside = true;
    for (size_t d=0; d<ndim; d++)
        inside = inside && pos[d]>=d_range[2*d] && pos[d]<=d_range[2*d+1];
    return inside;
}


/********************************************************
* Displace the mesh                                     *
********************************************************/
void Box::displaceMesh( const std::vector<double> &x )
{
    for (size_t d=0; d<d_range.size()/2; d++) {
        d_range[2*d+0] += x[d];
        d_range[2*d+1] += x[d];
    }
}


} // Geometry namespace
} // AMP namespace

