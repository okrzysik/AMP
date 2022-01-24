// This file stores the routines to generate the meshes for AMP::Mesh::Mesh
#include "AMP/TPLs.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshPoint.h"
#include "AMP/mesh/MultiMesh.h"
#include "AMP/mesh/structured/BoxMesh.h"
#include "AMP/mesh/triangle/TriangleHelpers.h"
#include "AMP/utils/Utilities.h"

#ifdef USE_TRILINOS_STKCLASSIC
//#include "AMP/mesh/STKmesh/STKMesh.h"
#endif
#ifdef USE_EXT_LIBMESH
    #include "AMP/mesh/libmesh/libmeshMesh.h"
#endif
#ifdef USE_EXT_MOAB
    #include "AMP/mesh/moab/moabMesh.h"
#endif


#include <cmath>


namespace AMP::Mesh {


std::map<std::string, AMP::Mesh::Mesh::generatorType> AMP::Mesh::Mesh::d_generators;


/********************************************************
 * Create a mesh from the input database                 *
 ********************************************************/
std::shared_ptr<AMP::Mesh::Mesh> Mesh::buildMesh( std::shared_ptr<const MeshParameters> params )
{
    auto db = params->d_db;
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
        auto suffix   = Utilities::getSuffix( filename );
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
#ifdef USE_EXT_LIBMESH
        mesh = std::make_shared<AMP::Mesh::libmeshMesh>( params );
#else
        AMP_ERROR( "AMP was compiled without support for libMesh" );
#endif
    } else if ( MeshType == "STKMesh" ) {
// The mesh is a stk mesh
#ifdef USE_TRILINOS_STKClassic
        // mesh = std::make_shared<AMP::Mesh::STKMesh>( params );
        AMP_ERROR( "AMP stk mesh interface is broken" );
#else
        AMP_ERROR( "AMP was compiled without support for STKMesh" );
#endif
    } else if ( MeshType == "moab" || MeshType == "MOAB" ) {
// The mesh is a MOAB mesh
#ifdef USE_EXT_MOAB
        mesh = std::make_shared<AMP::Mesh::moabMesh>( params );
#else
        AMP_ERROR( "AMP was compiled without support for MOAB" );
#endif
    } else {
        // Search for a mesh generator
        auto it = d_generators.find( MeshType );
        if ( it != d_generators.end() ) {
            mesh = it->second( params );
        } else {
            // Unknown mesh type
            AMP_ERROR( "Unknown mesh type (" + MeshType + ")" );
        }
    }
    mesh->setName( MeshName );
    return mesh;
}


/********************************************************
 * Estimate the mesh size                                *
 ********************************************************/
size_t Mesh::estimateMeshSize( std::shared_ptr<const MeshParameters> params )
{
    auto db = params->d_db;
    AMP_ASSERT( db != nullptr );
    size_t meshSize = 0;
    if ( db->keyExists( "NumberOfElements" ) ) {
        // User specified the number of elements, this should override everything
        meshSize = (size_t) db->getScalar<int>( "NumberOfElements" );
        // Adjust the number of elements by a weight if desired
        if ( db->keyExists( "Weight" ) ) {
            auto weight = db->getScalar<double>( "Weight" );
            meshSize    = (size_t) ceil( weight * ( (double) meshSize ) );
        }
        return meshSize;
    }
    // This is being called through the base class, call the appropriate function
    AMP_INSIST( db->keyExists( "MeshType" ), "MeshType must exist in input database" );
    auto MeshType = db->getString( "MeshType" );
    if ( db->keyExists( "NumberOfElements" ) ) {
        meshSize = db->getScalar<int>( "NumberOfElements" );
    } else if ( MeshType == "Multimesh" ) {
        // The mesh is a multimesh
        meshSize = AMP::Mesh::MultiMesh::estimateMeshSize( params );
    } else if ( MeshType == "AMP" ) {
        // The mesh is a AMP mesh
        auto filename = db->getWithDefault<std::string>( "FileName", "" );
        auto suffix   = Utilities::getSuffix( filename );
        if ( suffix == "stl" ) {
            // We are reading an stl file
            meshSize = AMP::Mesh::TriangleHelpers::readSTLHeader( filename );
        } else {
            meshSize = AMP::Mesh::BoxMesh::estimateMeshSize( params );
        }
    } else if ( MeshType == "TriangleGeometryMesh" ) {
        // We will build a triangle mesh from a geometry
        auto geom_db    = db->getDatabase( "Geometry" );
        auto geometry   = AMP::Geometry::Geometry::buildGeometry( geom_db );
        auto [lb, ub]   = geometry->box();
        auto resolution = db->getScalar<double>( "Resolution" );
        meshSize        = 1;
        for ( int d = 0; d < lb.ndim(); d++ )
            meshSize *= std::max<int64_t>( ( ub[d] - lb[d] ) / resolution, 1 );
    } else if ( MeshType == "libMesh" ) {
// The mesh is a libmesh mesh
#ifdef USE_EXT_LIBMESH
        meshSize = AMP::Mesh::libmeshMesh::estimateMeshSize( params );
#else
        AMP_ERROR( "AMP was compiled without support for libMesh" );
#endif
    } else if ( MeshType == "STKMesh" ) {
// The mesh is a stkMesh mesh
#ifdef USE_TRILINOS_STKCLASSIC
        // meshSize = AMP::Mesh::STKMesh::estimateMeshSize( params );
        AMP_ERROR( "AMP stk mesh interface is broken" );
#else
        AMP_ERROR( "AMP was compiled without support for STKMesh" );
#endif
    } else {
        // Unknown mesh type
        AMP_ERROR( "Unknown mesh type " + MeshType + " and NumberOfElements is not set" );
    }
    return meshSize;
}


/********************************************************
 * Estimate the maximum number of processors             *
 ********************************************************/
size_t Mesh::maxProcs( std::shared_ptr<const MeshParameters> params )
{
    auto db = params->d_db;
    AMP_ASSERT( db != nullptr );
    // Check if the user is specifying the maximum number of processors
    if ( db->keyExists( "maxProcs" ) )
        return db->getScalar<int64_t>( "maxProcs" );
    // This is being called through the base class, call the appropriate function
    AMP_INSIST( db->keyExists( "MeshType" ), "MeshType must exist in input database" );
    std::string MeshType = db->getString( "MeshType" );
    size_t maxSize       = 0;
    if ( db->keyExists( "maxProcs" ) ) {
        maxSize = db->getScalar<int>( "maxProcs" );
    } else if ( MeshType == std::string( "Multimesh" ) ) {
        // The mesh is a multimesh
        maxSize = AMP::Mesh::MultiMesh::maxProcs( params );
    } else if ( MeshType == std::string( "AMP" ) ) {
        // The mesh is a AMP mesh
        auto filename = db->getWithDefault<std::string>( "FileName", "" );
        auto suffix   = Utilities::getSuffix( filename );
        if ( suffix == "stl" ) {
            // We are reading an stl file
            maxSize = AMP::Mesh::TriangleHelpers::readSTLHeader( filename );
        } else {
            maxSize = AMP::Mesh::BoxMesh::maxProcs( params );
        }
    } else if ( MeshType == "TriangleGeometryMesh" ) {
        // We will build a triangle mesh from a geometry
        /*auto geom_db  = db->getDatabase( "Geometry" );
        auto geometry   = AMP::Geometry::Geometry::buildGeometry( geom_db );
        auto [lb, ub]   = geometry->box();
        auto resolution = db->getScalar<double>( "Resolution" );
        meshSize        = 1;
        for ( int d = 0; d < lb.ndim(); d++ )
            meshSize *= std::max<int64_t>( ( ub[d] - lb[d] ) / resolution, 1 ); */
        maxSize = 1;
    } else if ( MeshType == std::string( "libMesh" ) ) {
// The mesh is a libmesh mesh
#ifdef USE_EXT_LIBMESH
        maxSize = AMP::Mesh::libmeshMesh::maxProcs( params );
#else
        AMP_ERROR( "AMP was compiled without support for libMesh" );
#endif
    } else if ( MeshType == std::string( "STKMesh" ) ) {
// The mesh is a stkMesh mesh
#ifdef USE_TRILINOS_STKCLASSIC
        // maxSize = AMP::Mesh::STKMesh::maxProcs( params );
        AMP_ERROR( "AMP stk mesh interface is broken" );
#else
        AMP_ERROR( "AMP was compiled without support for STKMesh" );
#endif
    } else {
        // Unknown mesh type
        AMP_ERROR( "Unknown mesh type " + MeshType + " and maxProcs is not set" );
    }
    return maxSize;
}


} // namespace AMP::Mesh
