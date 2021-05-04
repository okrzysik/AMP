#include "AMP/ampmesh/Geometry.h"
#include "AMP/ampmesh/LogicalGeometry.h"
#include "AMP/ampmesh/MultiMesh.h"
#include "AMP/ampmesh/testHelpers/geometryTests.h"
#include "AMP/ampmesh/testHelpers/meshTests.h"
#include "AMP/utils/UnitTest.h"


namespace AMP {
namespace Mesh {


// This runs the geometry only tests
void meshTests::TestBasicGeometry( AMP::UnitTest &ut, AMP::Mesh::Mesh::const_shared_ptr mesh )
{
    // If we are dealing with a MultiMesh, check each mesh independently
    if ( std::dynamic_pointer_cast<const AMP::Mesh::MultiMesh>( mesh ) ) {
        auto multimesh = std::dynamic_pointer_cast<const AMP::Mesh::MultiMesh>( mesh );
        for ( const auto &mesh2 : multimesh->getMeshes() )
            TestBasicGeometry( ut, mesh2 );
    }
    // Get the geometry
    auto geom = mesh->getGeometry();
    if ( !geom )
        return;
    // Run basic geometry tests
    Geometry::testGeometry( *geom, ut );
}


// This tests loops over the boundary
void meshTests::TestInside( AMP::UnitTest &ut, AMP::Mesh::Mesh::const_shared_ptr mesh )
{
    // If we are dealing with a MultiMesh, check each mesh independently
    if ( std::dynamic_pointer_cast<const AMP::Mesh::MultiMesh>( mesh ) ) {
        auto multimesh = std::dynamic_pointer_cast<const AMP::Mesh::MultiMesh>( mesh );
        for ( const auto &mesh2 : multimesh->getMeshes() )
            TestInside( ut, mesh2 );
        return; // Eventually this should go away to test the multigeometry
    }
    // Get the geometry
    auto geom = mesh->getGeometry();
    if ( !geom )
        return;
    // Verify all elements in the mesh are inside the geometry
    bool pass                              = true;
    int gcw                                = mesh->getMaxGhostWidth();
    std::vector<AMP::Mesh::GeomType> types = { AMP::Mesh::GeomType::Vertex };
    if ( mesh->getDim() == static_cast<int>( mesh->getGeomType() ) )
        types.push_back( mesh->getGeomType() );
    for ( auto type : types ) {
        for ( const auto &elem : mesh->getIterator( type, gcw ) ) {
            auto p = elem.centroid();
            pass   = pass && geom->inside( p );
        }
    }
    if ( pass )
        ut.passes( "All mesh elements are inside geometry: " + mesh->getName() );
    else
        ut.failure( "Mesh elements are outside geometry: " + mesh->getName() );
}


// This tests checks physical-logical-physical transformations
void meshTests::TestPhysicalLogical( AMP::UnitTest &ut, AMP::Mesh::Mesh::const_shared_ptr mesh )
{
    bool pass = true;
    // If we are dealing with a MultiMesh, check each mesh independently
    if ( std::dynamic_pointer_cast<const AMP::Mesh::MultiMesh>( mesh ) ) {
        auto multimesh = std::dynamic_pointer_cast<const AMP::Mesh::MultiMesh>( mesh );
        for ( const auto &mesh2 : multimesh->getMeshes() )
            TestPhysicalLogical( ut, mesh2 );
        return;
    }
    // Get the geometry
    auto geom = std::dynamic_pointer_cast<AMP::Geometry::LogicalGeometry>( mesh->getGeometry() );
    if ( !geom )
        return;
    // Check the transformation of all points in the mesh
    auto type = AMP::Mesh::GeomType::Vertex;
    auto it   = mesh->getIterator( type, 0 );
    for ( size_t i = 0; i < it.size(); i++, ++it ) {
        auto p  = it->centroid();
        auto l  = geom->logical( p );
        auto p2 = geom->physical( l );
        bool t  = fabs( p.x() - p2.x() ) < 1e-7 && fabs( p.y() - p2.y() ) < 1e-7 &&
                 fabs( p.z() - p2.z() ) < 1e-7;
        pass = pass && t;
    }
    // Store the results
    if ( pass )
        ut.passes( "physical-logical-physical: " + mesh->getName() );
    else
        ut.failure( "physical-logical-physical: " + mesh->getName() );
}


} // namespace Mesh
} // namespace AMP
