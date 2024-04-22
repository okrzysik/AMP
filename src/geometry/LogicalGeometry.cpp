#include "AMP/geometry/LogicalGeometry.h"
#include "AMP/IO/HDF5.h"

#include <vector>


namespace AMP::Geometry {


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
