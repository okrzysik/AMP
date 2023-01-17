#include "AMP/mesh/MeshFactory.h"
#include "AMP/IO/FileSystem.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/mesh/MeshPoint.h"
#include "AMP/mesh/MultiMesh.h"
#include "AMP/mesh/structured/BoxMesh.h"
#include "AMP/mesh/triangle/TriangleHelpers.h"

#ifdef AMP_USE_TRILINOS_STKCLASSIC
//#include "AMP/mesh/STKmesh/STKMesh.h"
#endif
#ifdef AMP_USE_LIBMESH
    #include "AMP/mesh/libmesh/libmeshMesh.h"
#endif
#ifdef AMP_USE_MOAB
    #include "AMP/mesh/moab/moabMesh.h"
#endif


namespace AMP::Mesh {


// Create the operator
std::shared_ptr<Mesh> MeshFactory::create( std::shared_ptr<MeshParameters> params )
{
    auto db = params->getDatabase();
    AMP_ASSERT( db != nullptr );
    AMP_INSIST( db->keyExists( "MeshType" ), "MeshType must exist in input database" );
    AMP_INSIST( db->keyExists( "MeshName" ), "MeshName must exist in input database" );
    auto MeshType = db->getString( "MeshType" );
    auto MeshName = db->getString( "MeshName" );
    std::shared_ptr<AMP::Mesh::Mesh> mesh;
    if ( MeshType == "Multimesh" ) {
        // The mesh is a multimesh
        mesh = std::make_shared<AMP::Mesh::MultiMesh>( params );
    } else if ( MeshType == "AMP" ) {
        // The mesh is a AMP mesh
        auto filename = db->getWithDefault<std::string>( "FileName", "" );
        auto suffix   = IO::getSuffix( filename );
        if ( suffix == "stl" ) {
            // We are reading an stl file
            mesh = AMP::Mesh::TriangleHelpers::generateSTL( params );
        } else {
            mesh = AMP::Mesh::BoxMesh::generate( params );
        }
    } else if ( MeshType == "TriangleGeometryMesh" ) {
        // We will build a triangle mesh from a geometry
        auto geom_db   = db->getDatabase( "Geometry" );
        double dist[3] = { db->getWithDefault<double>( "x_offset", 0.0 ),
                           db->getWithDefault<double>( "y_offset", 0.0 ),
                           db->getWithDefault<double>( "z_offset", 0.0 ) };
        auto geom      = AMP::Geometry::Geometry::buildGeometry( geom_db );
        geom->displace( dist );
        auto res = db->getScalar<double>( "Resolution" );
        mesh     = AMP::Mesh::TriangleHelpers::generate( geom, params->getComm(), res );
    } else if ( MeshType == "libMesh" ) {
// The mesh is a libmesh mesh
#ifdef AMP_USE_LIBMESH
        mesh = std::make_shared<AMP::Mesh::libmeshMesh>( params );
#else
        AMP_ERROR( "AMP was compiled without support for libMesh" );
#endif
    } else if ( MeshType == "STKMesh" ) {
// The mesh is a stk mesh
#ifdef AMP_USE_TRILINOS_STKClassic
        // mesh = std::make_shared<AMP::Mesh::STKMesh>( params );
        AMP_ERROR( "AMP stk mesh interface is broken" );
#else
        AMP_ERROR( "AMP was compiled without support for STKMesh" );
#endif
    } else if ( MeshType == "moab" || MeshType == "MOAB" ) {
// The mesh is a MOAB mesh
#ifdef AMP_USE_MOAB
        mesh = std::make_shared<AMP::Mesh::moabMesh>( params );
#else
        AMP_ERROR( "AMP was compiled without support for MOAB" );
#endif
    } else {
        // Search for a mesh generator
        FactoryStrategy<Mesh, std::shared_ptr<MeshParameters>>::create( MeshType, params );
    }
    mesh->setName( MeshName );
    return mesh;
}


} // namespace AMP::Mesh
