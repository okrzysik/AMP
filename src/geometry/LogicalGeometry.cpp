#include "AMP/geometry/LogicalGeometry.h"

#include <vector>


namespace AMP::Geometry {


/********************************************************
 * LogicalGeometry                                       *
 ********************************************************/
std::vector<bool> LogicalGeometry::getPeriodicDim() const
{
    return std::vector<bool>( d_isPeriodic.data(), d_isPeriodic.data() + d_logicalDim );
}
std::vector<int> LogicalGeometry::getLogicalSurfaceIds() const
{
    return std::vector<int>( d_ids.data(), d_ids.data() + 2 * d_logicalDim );
}

} // namespace AMP::Geometry
