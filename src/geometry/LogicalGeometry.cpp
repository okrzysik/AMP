#include "AMP/geometry/LogicalGeometry.h"
#include "AMP/IO/HDF5.h"

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


/********************************************************
 * Write/Read restart data                               *
 ********************************************************/
void LogicalGeometry::writeRestart( int64_t fid ) const
{
    Geometry::writeRestart( fid );
    writeHDF5( fid, "logicalDim", d_logicalDim );
    writeHDF5( fid, "isPeriodic", d_isPeriodic );
    writeHDF5( fid, "ids", d_ids );
}
LogicalGeometry::LogicalGeometry( int64_t fid ) : Geometry( fid )
{
    readHDF5( fid, "logicalDim", d_logicalDim );
    readHDF5( fid, "isPeriodic", d_isPeriodic );
    readHDF5( fid, "ids", d_ids );
}


} // namespace AMP::Geometry
