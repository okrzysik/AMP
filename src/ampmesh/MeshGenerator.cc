// This file stores the routines to generate the meshes for AMP::Mesh::Mesh
#include "AMP/ampmesh/Mesh.h"
#include "AMP/utils/Utilities.h"

#include "AMP/ampmesh/MultiMesh.h"
#include "AMP/ampmesh/structured/BoxMesh.h"
#include "AMP/ampmesh/triangle/TriangleMesh.h"
#ifdef USE_TRILINOS_STKCLASSIC
//#include "AMP/ampmesh/STKmesh/STKMesh.h"
#endif
#ifdef USE_EXT_LIBMESH
#include "AMP/ampmesh/libmesh/libMesh.h"
#endif
#ifdef USE_EXT_MOAB
#include "AMP/ampmesh/moab/moabMesh.h"
#endif


#include <cmath>


namespace AMP {
namespace Mesh {


std::map<std::string, AMP::Mesh::Mesh::generatorType> AMP::Mesh::Mesh::d_generators;


/********************************************************
 * Create a mesh from the input database                 *
 ********************************************************/
AMP::shared_ptr<AMP::Mesh::Mesh> Mesh::buildMesh( const MeshParameters::shared_ptr &params )
{
    auto database = params->d_db;
    AMP_ASSERT( database != nullptr );
    AMP_INSIST( database->keyExists( "MeshType" ), "MeshType must exist in input database" );
    AMP_INSIST( database->keyExists( "MeshName" ), "MeshName must exist in input database" );
    std::string MeshType = database->getString( "MeshType" );
    std::string MeshName = database->getString( "MeshName" );
    AMP::shared_ptr<AMP::Mesh::Mesh> mesh;
    if ( MeshType == std::string( "Multimesh" ) ) {
        // The mesh is a multimesh
        mesh = AMP::make_shared<AMP::Mesh::MultiMesh>( params );
    } else if ( MeshType == std::string( "AMP" ) ) {
        // The mesh is a AMP mesh
        auto filename = database->getWithDefault<std::string>( "FileName", "" );
        if ( filename.substr( std::max<int>( filename.length(), 4 ) - 4 ) == ".stl" ) {
            // We are reading an stl file
            mesh = AMP::Mesh::TriangleMesh<2, 3>::generate( params );
        } else {
            mesh = AMP::Mesh::BoxMesh::generate( params );
        }
    } else if ( MeshType == std::string( "libMesh" ) ) {
// The mesh is a libmesh mesh
#ifdef USE_EXT_LIBMESH
        mesh = AMP::make_shared<AMP::Mesh::libMesh>( params );
#else
        AMP_ERROR( "AMP was compiled without support for libMesh" );
#endif
    } else if ( MeshType == std::string( "STKMesh" ) ) {
// The mesh is a stk mesh
#ifdef USE_TRILINOS_STKClassic
        // mesh = AMP::make_shared<AMP::Mesh::STKMesh>( params );
        AMP_ERROR( "AMP stk mesh interface is broken" );
#else
        AMP_ERROR( "AMP was compiled without support for STKMesh" );
#endif
    } else if ( MeshType == std::string( "moab" ) || MeshType == std::string( "MOAB" ) ) {
// The mesh is a MOAB mesh
#ifdef USE_EXT_MOAB
        mesh = AMP::make_shared<AMP::Mesh::moabMesh>( params );
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
            AMP_ERROR( std::string( "Unknown mesh type (" ) + MeshType + std::string( ")" ) );
        }
    }
    mesh->setName( MeshName );
    return mesh;
}


/********************************************************
 * Estimate the mesh size                                *
 ********************************************************/
size_t Mesh::estimateMeshSize( const MeshParameters::shared_ptr &params )
{
    auto database = params->d_db;
    AMP_ASSERT( database != nullptr );
    size_t meshSize = 0;
    if ( database->keyExists( "NumberOfElements" ) ) {
        // User specified the number of elements, this should override everything
        meshSize = (size_t) database->getScalar<int>( "NumberOfElements" );
        // Adjust the number of elements by a weight if desired
        if ( database->keyExists( "Weight" ) ) {
            double weight = database->getScalar<double>( "Weight" );
            meshSize      = (size_t) ceil( weight * ( (double) meshSize ) );
        }
        return meshSize;
    }
    // This is being called through the base class, call the appropriate function
    AMP_INSIST( database->keyExists( "MeshType" ), "MeshType must exist in input database" );
    std::string MeshType = database->getString( "MeshType" );
    if ( MeshType == std::string( "Multimesh" ) ) {
        // The mesh is a multimesh
        meshSize = AMP::Mesh::MultiMesh::estimateMeshSize( params );
    } else if ( MeshType == std::string( "AMP" ) ) {
        // The mesh is a AMP mesh
        auto filename = database->getWithDefault<std::string>( "FileName", "" );
        if ( filename.substr( std::max<int>( filename.length(), 4 ) - 4 ) == ".stl" ) {
            // We are reading an stl file
            meshSize = AMP::Mesh::TriangleMesh<2, 3>::estimateMeshSize( params );
        } else {
            meshSize = AMP::Mesh::BoxMesh::estimateMeshSize( params );
        }
    } else if ( MeshType == std::string( "libMesh" ) ) {
// The mesh is a libmesh mesh
#ifdef USE_EXT_LIBMESH
        meshSize = AMP::Mesh::libMesh::estimateMeshSize( params );
#else
        AMP_ERROR( "AMP was compiled without support for libMesh" );
#endif
    } else if ( MeshType == std::string( "STKMesh" ) ) {
// The mesh is a stkMesh mesh
#ifdef USE_TRILINOS_STKCLASSIC
        // meshSize = AMP::Mesh::STKMesh::estimateMeshSize( params );
        AMP_ERROR( "AMP stk mesh interface is broken" );
#else
        AMP_ERROR( "AMP was compiled without support for STKMesh" );
#endif
    } else if ( database->keyExists( "NumberOfElements" ) ) {
        int NumberOfElements = database->getScalar<int>( "NumberOfElements" );
        meshSize             = NumberOfElements;
    } else {
        // Unknown mesh type
        AMP_ERROR( "Unknown mesh type and NumberOfElements does not exist in database" );
    }
    return meshSize;
}


/********************************************************
 * Estimate the maximum number of processors             *
 ********************************************************/
size_t Mesh::maxProcs( const MeshParameters::shared_ptr &params )
{
    auto database = params->d_db;
    AMP_ASSERT( database != nullptr );
    // This is being called through the base class, call the appropriate function
    AMP_INSIST( database->keyExists( "MeshType" ), "MeshType must exist in input database" );
    std::string MeshType = database->getString( "MeshType" );
    size_t maxSize       = 0;
    if ( MeshType == std::string( "Multimesh" ) ) {
        // The mesh is a multimesh
        maxSize = AMP::Mesh::MultiMesh::maxProcs( params );
    } else if ( MeshType == std::string( "AMP" ) ) {
        // The mesh is a AMP mesh
        maxSize = AMP::Mesh::BoxMesh::maxProcs( params );
    } else if ( MeshType == std::string( "libMesh" ) ) {
// The mesh is a libmesh mesh
#ifdef USE_EXT_LIBMESH
        maxSize = AMP::Mesh::libMesh::maxProcs( params );
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
        AMP_ERROR( "Unknown mesh type and NumberOfElements does not exist in database" );
    }
    return maxSize;
}


} // namespace Mesh
} // namespace AMP
