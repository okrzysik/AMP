// This file contains classes for generating meshes that are based in libmesh
#ifndef included_AMP_Unit_test_Libmesh_Generators_h
#define included_AMP_Unit_test_Libmesh_Generators_h

#include "AMP/ampmesh/libmesh/initializeLibMesh.h"
#include "AMP/ampmesh/libmesh/libmeshMesh.h"
#include "meshGenerators.h"

// LibMesh include
#include "libmesh/cell_hex8.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"


namespace AMP {
namespace unit_test {


// Class to create a cube in Libmesh
template<int SIZE>
class LibMeshCubeGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() override
    {
        // Create the parameter object
        auto database = std::make_shared<AMP::Database>( "Mesh" );
        database->putScalar<int>( "dim", 3 );
        database->putScalar<std::string>( "MeshName", "cube_mesh" );
        database->putScalar<std::string>( "Generator", "cube" );
        database->putVector<int>( "size", std::vector<int>( 3, SIZE ) );
        database->putVector<double>( "xmin", std::vector<double>( 3, -1.0 ) );
        database->putVector<double>( "xmax", std::vector<double>( 3, 1.0 ) );
        auto params = std::make_shared<AMP::Mesh::MeshParameters>( database );
        params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        // Create a libMesh mesh
        mesh = std::make_shared<AMP::Mesh::libmeshMesh>( params );
    }

    static std::string name() { return "LibMeshCubeGenerator"; }
};


// Class to read in a default exodus file
template<int FILE = 1>
class ExodusReaderGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() override
    {
        // Create the parameter object
        auto database = std::make_shared<AMP::Database>( "Mesh" );
        database->putScalar( "dim", 3 );
        database->putScalar<std::string>( "MeshName", "exodus reader mesh" );
        if ( FILE == 1 ) {
            database->putScalar<std::string>( "FileName", "clad_1x_1pellet.e" );
        } else if ( FILE == 2 ) {
            database->putScalar<std::string>( "FileName", "multiElementMesh.e" );
        } else if ( FILE == 3 ) {
            database->putScalar<std::string>( "FileName", "pellet_1x.e" );
        } else {
            AMP_ERROR( "Bad file for generator" );
        }
        auto params = std::make_shared<AMP::Mesh::MeshParameters>( database );
        params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        // Create a libMesh mesh
        mesh = std::make_shared<AMP::Mesh::libmeshMesh>( params );
    }

    static std::string name() { return "ExodusReaderGenerator"; }
};


// MulitMesh generator
class MultiMeshGenerator : public MeshGenerator
{
public:
    virtual void build_mesh() override
    {
        int N_meshes = 4;
        // Create the multimesh database
        auto meshDatabase = std::make_shared<AMP::Database>( "Mesh" );
        meshDatabase->putScalar<std::string>( "MeshName", "PelletMeshes" );
        meshDatabase->putScalar<std::string>( "MeshType", "Multimesh" );
        meshDatabase->putScalar<std::string>( "MeshDatabasePrefix", "Mesh_" );
        meshDatabase->putScalar<std::string>( "MeshArrayDatabasePrefix", "MeshArray_" );
        // Create the mesh array database
        auto meshArrayDatabase = meshDatabase->putDatabase( "MeshArray_1" );
        meshArrayDatabase->putScalar<int>( "N", N_meshes );
        meshArrayDatabase->putScalar<std::string>( "iterator", "%i" );
        std::vector<int> indexArray( N_meshes );
        for ( int i = 0; i < N_meshes; i++ )
            indexArray[i] = i + 1;
        meshArrayDatabase->putVector<int>( "indicies", indexArray );
        meshArrayDatabase->putScalar<std::string>( "MeshName", "pellet_%i" );
#ifdef USE_EXT_LIBMESH
        meshArrayDatabase->putScalar<std::string>( "FileName", "pellet_lo_res.e" );
        meshArrayDatabase->putScalar<std::string>( "MeshType", "libMesh" );
#else
        std::vector<int> size( 3, 10 );
        std::vector<double> range( 6, 0.0 );
        range[1] = 0.005;
        range[3] = 0.005;
        range[5] = 0.005;
        // Create a generic MeshParameters object
        meshArrayDatabase->putScalar<std::string>( "MeshType", "AMP" );
        meshArrayDatabase->putScalar<std::string>( "Generator", "cube" );
        meshArrayDatabase->putVector<int>( "Size", size );
        meshArrayDatabase->putVector<double>( "Range", range );
#endif
        meshArrayDatabase->putScalar<double>( "dim", 3 );
        meshArrayDatabase->putScalar<double>( "x_offset", 0.0 );
        meshArrayDatabase->putScalar<double>( "y_offset", 0.0 );
        std::vector<double> offsetArray( N_meshes );
        for ( int i = 0; i < N_meshes; i++ )
            offsetArray[i] = ( (double) i ) * 0.0105;
        meshArrayDatabase->putVector<double>( "z_offset", offsetArray );
        // Create the parameter object
        auto params = std::make_shared<AMP::Mesh::MeshParameters>( meshDatabase );
        params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        // Create the mesh
        mesh = AMP::Mesh::Mesh::buildMesh( params );
    }
    static std::string name() { return "MultiMeshGenerator"; }
};


// libMeshThreeElement generator
class libMeshThreeElementGenerator : public MeshGenerator
{
public:
    static std::string name() { return "libMeshThreeElementGenerator"; }

    static std::vector<unsigned int> getBndDofIndices()
    {
        std::vector<unsigned int> bndDofIndices( 4 );
        bndDofIndices[0] = 0;
        bndDofIndices[1] = 3;
        bndDofIndices[2] = 4;
        bndDofIndices[3] = 7;
        return bndDofIndices;
    }

    static std::vector<std::vector<unsigned int>> getElemNodeMap()
    {
        std::vector<std::vector<unsigned int>> elemNodeMap( 3 );

        elemNodeMap[0].resize( 8 );
        elemNodeMap[1].resize( 8 );
        elemNodeMap[2].resize( 8 );

        elemNodeMap[0][0] = 0;
        elemNodeMap[0][1] = 1;
        elemNodeMap[0][2] = 2;
        elemNodeMap[0][3] = 3;
        elemNodeMap[0][4] = 4;
        elemNodeMap[0][5] = 5;
        elemNodeMap[0][6] = 6;
        elemNodeMap[0][7] = 7;

        elemNodeMap[1][0] = 1;
        elemNodeMap[1][1] = 8;
        elemNodeMap[1][2] = 9;
        elemNodeMap[1][3] = 2;
        elemNodeMap[1][4] = 5;
        elemNodeMap[1][5] = 10;
        elemNodeMap[1][6] = 11;
        elemNodeMap[1][7] = 6;

        elemNodeMap[2][0] = 2;
        elemNodeMap[2][1] = 9;
        elemNodeMap[2][2] = 12;
        elemNodeMap[2][3] = 13;
        elemNodeMap[2][4] = 6;
        elemNodeMap[2][5] = 11;
        elemNodeMap[2][6] = 14;
        elemNodeMap[2][7] = 15;

        return elemNodeMap;
    }

    virtual void build_mesh() override
    {

        // Initialize libmesh
        AMP::AMP_MPI comm( AMP_COMM_SELF );
        libmeshInit = std::make_shared<AMP::Mesh::initializeLibMesh>( comm );
        libMeshComm = std::make_shared<libMesh::Parallel::Communicator>( comm.getCommunicator() );

        const unsigned int mesh_dim  = 3;
        const unsigned int num_elem  = 3;
        const unsigned int num_nodes = 16;

        auto local_mesh = std::make_shared<libMesh::Mesh>( *libMeshComm, mesh_dim );
        local_mesh->reserve_elem( num_elem );
        local_mesh->reserve_nodes( num_nodes );

        local_mesh->add_point(::libMesh::Point( 0.0, 0.0, 0.0 ), 0 );
        local_mesh->add_point(::libMesh::Point( 0.5, 0.0, 0.0 ), 1 );
        local_mesh->add_point(::libMesh::Point( 0.5, 0.5, 0.0 ), 2 );
        local_mesh->add_point(::libMesh::Point( 0.0, 0.5, 0.0 ), 3 );
        local_mesh->add_point(::libMesh::Point( 0.0, 0.0, 0.5 ), 4 );
        local_mesh->add_point(::libMesh::Point( 0.5, 0.0, 0.5 ), 5 );
        local_mesh->add_point(::libMesh::Point( 0.5, 0.5, 0.5 ), 6 );
        local_mesh->add_point(::libMesh::Point( 0.0, 0.5, 0.5 ), 7 );
        local_mesh->add_point(::libMesh::Point( 1.0, 0.0, 0.0 ), 8 );
        local_mesh->add_point(::libMesh::Point( 1.0, 0.5, 0.0 ), 9 );
        local_mesh->add_point(::libMesh::Point( 1.0, 0.0, 0.5 ), 10 );
        local_mesh->add_point(::libMesh::Point( 1.0, 0.5, 0.5 ), 11 );
        local_mesh->add_point(::libMesh::Point( 1.0, 1.0, 0.0 ), 12 );
        local_mesh->add_point(::libMesh::Point( 0.5, 1.0, 0.0 ), 13 );
        local_mesh->add_point(::libMesh::Point( 1.0, 1.0, 0.5 ), 14 );
        local_mesh->add_point(::libMesh::Point( 0.5, 1.0, 0.5 ), 15 );

        auto elemNodeMap = getElemNodeMap();
        for ( size_t i = 0; i < elemNodeMap.size(); i++ ) {
            libMesh::Elem *elem = local_mesh->add_elem( new ::libMesh::Hex8 );
            for ( int j = 0; j < 8; j++ ) {
                elem->set_node( j ) = local_mesh->node_ptr( elemNodeMap[i][j] );
            }
        }

        const short int boundaryId = 1;
        for ( auto &bndDofIndice : getBndDofIndices() ) {
            local_mesh->boundary_info->add_node( local_mesh->node_ptr( bndDofIndice ), boundaryId );
        }

        local_mesh->prepare_for_use( true );
        mesh = std::make_shared<AMP::Mesh::libmeshMesh>( local_mesh, "3 Element" );
    }

    virtual ~libMeshThreeElementGenerator()
    {
        mesh.reset();
        libmeshInit.reset();
    }

protected:
    std::shared_ptr<libMesh::Parallel::Communicator> libMeshComm;
    std::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit;
};


} // namespace unit_test
} // namespace AMP


#endif
