#include "AMP/geometry/Geometry.h"
#include "AMP/geometry/LogicalGeometry.h"
#include "AMP/geometry/shapes/RegularPolygon.h"
#include "AMP/geometry/testHelpers/geometryTests.h"
#include "AMP/mesh/MultiMesh.h"
#include "AMP/mesh/testHelpers/meshTests.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"


namespace AMP::Mesh {


// This runs the geometry only tests
void meshTests::TestBasicGeometry( AMP::UnitTest &ut, std::shared_ptr<const AMP::Mesh::Mesh> mesh )
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
void meshTests::TestInside( AMP::UnitTest &ut, std::shared_ptr<const AMP::Mesh::Mesh> mesh )
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
            if ( !geom->inside( p ) ) {
                p       = elem.centroid();
                bool in = geom->inside( p );
                NULL_USE( in );
                pass = false;
            }
        }
    }
    if ( pass )
        ut.passes( "All mesh elements are inside geometry: " + mesh->getName() );
    else
        ut.failure( "Mesh elements are outside geometry: " + mesh->getName() );
}


// This tests checks physical-logical-physical transformations
void meshTests::TestPhysicalLogical( AMP::UnitTest &ut,
                                     std::shared_ptr<const AMP::Mesh::Mesh> mesh )
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


// This tests that the normal in the geometry and the mesh agree (within discretization error)
void meshTests::TestNormalGeometry( AMP::UnitTest &ut, std::shared_ptr<const AMP::Mesh::Mesh> mesh )
{
    if ( mesh->getDim() <= 1 )
        return; // Normals are not defined for 1D
    // If we are dealing with a MultiMesh, check each mesh independently
    if ( std::dynamic_pointer_cast<const AMP::Mesh::MultiMesh>( mesh ) ) {
        auto multimesh = std::dynamic_pointer_cast<const AMP::Mesh::MultiMesh>( mesh );
        for ( const auto &mesh2 : multimesh->getMeshes() )
            TestNormalGeometry( ut, mesh2 );
        return;
    }
    // Get the geometry
    auto geom = mesh->getGeometry();
    if ( !geom )
        return;
    // Loop over the surface
    double error = 0.0;
    auto type    = static_cast<AMP::Mesh::GeomType>( mesh->getDim() - 1 );
    for ( auto &elem : mesh->getSurfaceIterator( type, 0 ) ) {
        // Get the normal from the geometry
        auto a  = elem.centroid();
        auto n1 = geom->surfaceNorm( a );
        // Get the normal from the element
        auto n2 = elem.norm();
        if ( dot( n1, n2 ) < 0 )
            n2 = -n2; // We do not always agree on direction (need to fix this)
        // Check the error
        auto err = abs( n1 - n2 );
        error    = std::max( error, err );
    }
    if ( error < 1e-6 ) {
        ut.passes( "mesh normal matches geom normal: " + mesh->getName() );
    } else {
        auto msg = AMP::Utilities::stringf( "%s (%f)", mesh->getName().data(), error );
        if ( error < 0.1 ) {
            ut.expected_failure( "mesh normal approximately matches geom normal: " + msg );
        } else if ( std::dynamic_pointer_cast<AMP::Geometry::RegularPolygon>( geom ) &&
                    error < 0.6 ) {
            // RegularPolygon has a larger error because the verticies of the polygon do
            //    not align with the verticies of the mesh
            ut.expected_failure( "mesh normal approximately matches geom normal: " + msg );
        } else {
            ut.failure( "mesh normal does not match geom normal: " + msg );
        }
    }
}


} // namespace AMP::Mesh
