#include "AMP/ampmesh/Geometry.h"
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshGeometry.h"
#include "AMP/ampmesh/MultiGeometry.h"
#include "AMP/ampmesh/MultiMesh.h"
#include "AMP/ampmesh/shapes/Box.h"
#include "AMP/ampmesh/shapes/Circle.h"
#include "AMP/ampmesh/shapes/CircleFrustum.h"
#include "AMP/ampmesh/shapes/Cylinder.h"
#include "AMP/ampmesh/shapes/Shell.h"
#include "AMP/ampmesh/shapes/Sphere.h"
#include "AMP/ampmesh/shapes/SphereSurface.h"
#include "AMP/ampmesh/shapes/SquareFrustum.h"
#include "AMP/ampmesh/shapes/Tube.h"
#include "AMP/utils/Database.h"


namespace AMP {
namespace Geometry {


/********************************************************
 * Create the geometry object                            *
 ********************************************************/
std::shared_ptr<AMP::Geometry::Geometry>
Geometry::buildGeometry( std::shared_ptr<AMP::Database> db )
{
    auto generator = db->getString( "Generator" );
    std::shared_ptr<AMP::Geometry::Geometry> geom;
    if ( generator.compare( "cube" ) == 0 ) {
        int dim = db->getScalar<int>( "dim" );
        if ( db->keyExists( "Range" ) ) {
            if ( dim == 1 ) {
                geom.reset( new Box<1>( db ) );
            } else if ( dim == 2 ) {
                geom.reset( new Box<2>( db ) );
            } else if ( dim == 3 ) {
                geom.reset( new Box<3>( db ) );
            } else {
                AMP_ERROR( "Physical Dimensions > 3 are not supported yet" );
            }
        } else if ( db->keyExists( "x_grid" ) ) {
            if ( dim == 1 ) {
                geom.reset( new Grid<1>( db ) );
            } else if ( dim == 2 ) {
                geom.reset( new Grid<2>( db ) );
            } else if ( dim == 3 ) {
                geom.reset( new Grid<3>( db ) );
            } else {
                AMP_ERROR( "Physical Dimensions > 3 are not supported yet" );
            }
        }
    } else if ( generator.compare( "tube" ) == 0 ) {
        geom.reset( new Tube( db ) );
    } else if ( generator.compare( "circle" ) == 0 ) {
        geom.reset( new Circle( db ) );
    } else if ( generator.compare( "cylinder" ) == 0 ) {
        geom.reset( new Cylinder( db ) );
    } else if ( generator.compare( "shell" ) == 0 ) {
        geom.reset( new Shell( db ) );
    } else if ( generator.compare( "sphere" ) == 0 ) {
        geom.reset( new Sphere( db ) );
    } else if ( generator.compare( "sphere_surface" ) == 0 ) {
        geom.reset( new SphereSurface( db ) );
    } else if ( generator.compare( "square_frustrum" ) == 0 ) {
        geom.reset( new SquareFrustum( db ) );
    } else if ( generator.compare( "circle_frustrum" ) == 0 ) {
        geom.reset( new CircleFrustum( db ) );
    } else if ( generator.compare( "Mesh" ) == 0 ) {
        // Generate a mesh geometry
        auto mesh_db = db->getDatabase( "Mesh" );
        auto params  = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
        params->setComm( AMP_COMM_SELF );
        auto mesh      = AMP::Mesh::Mesh::buildMesh( params );
        auto multimesh = std::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( mesh );
        if ( multimesh ) {
            std::vector<Geometry::shared_ptr> geoms;
            for ( const auto &mesh : multimesh->getMeshes() )
                geoms.push_back( std::make_shared<MeshGeometry>( mesh ) );
            geom = std::make_shared<MultiGeometry>( geoms );
        } else {
            geom = std::make_shared<MeshGeometry>( mesh );
        }
    } else {
        AMP_ERROR( "Unknown generator" );
    }
    return geom;
}


} // namespace Geometry
} // namespace AMP
