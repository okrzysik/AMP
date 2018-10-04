#include "AMP/ampmesh/Geometry.h"
#include "AMP/ampmesh/testHelpers/meshTests.h"
#include "AMP/utils/UnitTest.h"


namespace AMP {
namespace Mesh {


// This tests loops over the boundary
void meshTests::TestInside( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh )
{
    auto geom = mesh->getGeometry();
    // Verify all elements in the mesh are inside the geometry
    bool pass = true;
    int gcw   = mesh->getMaxGhostWidth();
    for ( int type2 = 0; type2 <= (int) mesh->getGeomType(); type2++ ) {
        auto type     = (AMP::Mesh::GeomType) type2;
        auto iterator = mesh->getIterator( type, gcw );
        for ( const auto &elem : mesh->getIterator( type, gcw ) ) {
            auto p = elem.centroid();
            pass   = pass && geom->inside( p );
        }
    }
    if ( pass )
        ut->passes( "All mesh elements are inside geometry" );
    else
        ut->failure( "Mesh elements are outside geometry" );
}


} // namespace Mesh
} // namespace AMP
