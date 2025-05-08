#include "AMP/mesh/libmesh/ReadTestMesh.h"
#include "AMP/IO/FileSystem.h"
#include "AMP/mesh/testHelpers/meshWriters.h"
#include "AMP/utils/Utilities.h"

// LibMesh include
DISABLE_WARNINGS
#include "libmesh/libmesh_config.h"
#undef LIBMESH_ENABLE_REFERENCE_COUNTING
#include "libmesh/boundary_info.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/elem.h"
#include "libmesh/mesh.h"
ENABLE_WARNINGS

#include <cstdio>

namespace AMP {

void readBinaryTestMesh( const std::string &mesh_file, std::shared_ptr<libMesh::Mesh> mesh )
{

    FILE *fp = fopen( mesh_file.c_str(), "rb" );

    int num_nodes;
    size_t n;
    n = fread( &num_nodes, sizeof( int ), 1, fp );
    AMP_INSIST( ( n == 1 ), "Error while reading the file" );

    mesh->reserve_nodes( num_nodes );

    std::vector<double> points( 3 * num_nodes );

    n = fread( &( points[0] ), sizeof( double ), ( 3 * num_nodes ), fp );
    AMP_INSIST( ( n == size_t( 3 * num_nodes ) ), "Error while reading the file" );

    for ( int i = 0; i < num_nodes; i++ ) {
        mesh->add_point(
            libMesh::Point( points[( 3 * i ) + 0], points[( 3 * i ) + 1], points[( 3 * i ) + 2] ),
            i );
    }

    points.clear();

    int num_elem;
    n = fread( &num_elem, sizeof( int ), 1, fp );
    AMP_INSIST( ( n == 1 ), "Error while reading the file" );

    mesh->reserve_elem( num_elem );

    std::vector<int> elemNodeMap( 8 * num_elem );

    n = fread( &( elemNodeMap[0] ), sizeof( int ), ( 8 * num_elem ), fp );
    AMP_INSIST( ( n == size_t( 8 * num_elem ) ), "Error while reading the file" );

    for ( int i = 0; i < num_elem; i++ ) {
        libMesh::Elem *newElem = new libMesh::Hex8;
        newElem->set_id( i );
        libMesh::Elem *elem = mesh->add_elem( newElem );
        for ( int j = 0; j < 8; j++ ) {
            elem->set_node( j ) = mesh->node_ptr( elemNodeMap[( 8 * i ) + j] );
        }
    }

    elemNodeMap.clear();

    int num_boundaryNodeIds;
    n = fread( &num_boundaryNodeIds, sizeof( int ), 1, fp );
    AMP_INSIST( ( n == 1 ), "Error while reading the file" );

    for ( int bid = 1; bid <= num_boundaryNodeIds; bid++ ) {
        int bndDofSize;
        n = fread( &bndDofSize, sizeof( int ), 1, fp );
        AMP_INSIST( ( n == 1 ), "Error while reading the file" );
        if ( bndDofSize > 0 ) {
            int idx;
            for ( int i = 0; i < bndDofSize; i++ ) {
                n = fread( &idx, sizeof( int ), 1, fp );
                AMP_INSIST( ( n == 1 ), "Error while reading the file" );
                mesh->boundary_info->add_node( mesh->node_ptr( idx ), bid );
            }
        }
    }

    int num_boundarySideIds;
    n = fread( &num_boundarySideIds, sizeof( int ), 1, fp );
    AMP_INSIST( ( n == 1 ), "Error while reading the file" );

    for ( int bid = 1; bid <= num_boundarySideIds; bid++ ) {
        int bndDofSize;
        n = fread( &bndDofSize, sizeof( int ), 1, fp );
        AMP_INSIST( ( n == 1 ), "Error while reading the file" );
        if ( bndDofSize > 0 ) {
            int idxE;
            int idxS;
            for ( int i = 0; i < bndDofSize; i++ ) {
                n = fread( &idxE, sizeof( int ), 1, fp );
                AMP_INSIST( ( n == 1 ), "Error while reading the file" );
                n = fread( &idxS, sizeof( int ), 1, fp );
                AMP_INSIST( ( n == 1 ), "Error while reading the file" );
                mesh->boundary_info->add_side( mesh->elem_ptr( idxE ), idxS, bid );
            }
        }
    }

    fclose( fp );
}

void readTestMesh( const std::string &mesh_file, std::shared_ptr<libMesh::Mesh> mesh )
{
    [[maybe_unused]] auto tmp = AMP::Utilities::randomString();
    if ( mesh_file == "distortedElementMesh" ) {
        AMP::Mesh::MeshWriters::writeDistortedElement( tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "cookMesh0" ) {
        AMP::Mesh::MeshWriters::writeCookMesh( 9, 2, 9, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "cookMesh1" ) {
        AMP::Mesh::MeshWriters::writeCookMesh( 17, 3, 17, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "cookMesh2" ) {
        AMP::Mesh::MeshWriters::writeCookMesh( 33, 5, 33, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "cookMesh3" ) {
        AMP::Mesh::MeshWriters::writeCookMesh( 65, 9, 65, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "cookMesh4" ) {
        AMP::Mesh::MeshWriters::writeCookMesh( 129, 17, 129, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "regPlateWithHole1" ) {
        AMP::Mesh::MeshWriters::writePlateWithHole( 5, 15, 15, 15, 5.0, 10.0, 1.0, 1.0, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "regPlateWithHole2" ) {
        AMP::Mesh::MeshWriters::writePlateWithHole( 10, 30, 30, 30, 5.0, 10.0, 1.0, 1.0, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "mesh7elem-1" ) {
        AMP::Mesh::MeshWriters::write7elementMesh( 2, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "mesh7elem-2" ) {
        AMP::Mesh::MeshWriters::write7elementMesh( 8, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "boxMesh-1" ) {
        AMP::Mesh::MeshWriters::writeAMGMesh( 2, 2, 2, 10, 10, 10, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "boxMesh-2" ) {
        AMP::Mesh::MeshWriters::writeAMGMesh( 3, 3, 3, 10, 10, 10, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "boxMesh-3" ) {
        AMP::Mesh::MeshWriters::writeAMGMesh( 5, 5, 5, 10, 10, 10, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "boxMesh-4" ) {
        AMP::Mesh::MeshWriters::writeAMGMesh( 9, 9, 9, 10, 10, 10, tmp );
        readTestMesh( tmp, mesh );
    } else if ( mesh_file == "boxMesh-5" ) {
        AMP::Mesh::MeshWriters::writeAMGMesh( 17, 17, 17, 10, 10, 10, tmp );
        readTestMesh( tmp, mesh );
    } else {
        auto db = AMP::Database::parseInputFile( mesh_file );
        readTestMesh( db, mesh );
    }
    AMP::IO::deleteFile( tmp );
}

void readTestMesh( std::shared_ptr<AMP::Database> db, std::shared_ptr<libMesh::Mesh> mesh )
{
    if ( db->keyExists( "Mesh" ) )
        db = db->getDatabase( "Mesh" );

    int num_elem            = db->getScalar<int>( "NumberOfElements" );
    int num_nodes           = db->getScalar<int>( "NumberOfNodes" );
    int num_boundaryNodeIds = db->getScalar<int>( "NumberOfBoundaryNodeIds" );
    int num_boundarySideIds = db->getScalar<int>( "NumberOfBoundarySideIds" );

    mesh->reserve_elem( num_elem );
    mesh->reserve_nodes( num_nodes );

    for ( int i = 0; i < num_nodes; i++ ) {
        char key[100];
        snprintf( key, sizeof key, "Point%d", i );
        auto point = db->getVector<double>( key );
        mesh->add_point( libMesh::Point( point[0], point[1], point[2] ), i );
    } // end for i

    std::vector<std::vector<int>> elemNodeMap;

    for ( int i = 0; i < num_elem; i++ ) {
        char key[100];
        snprintf( key, sizeof key, "Elem%d", i );
        elemNodeMap.push_back( db->getVector<int>( key ) );
    } // end for i

    for ( int i = 0; i < num_elem; i++ ) {
        libMesh::Elem *elem = mesh->add_elem( new libMesh::Hex8 );
        for ( int j = 0; j < 8; j++ ) {
            elem->set_node( j ) = mesh->node_ptr( elemNodeMap[i][j] );
        } // end for j
    }     // end for i

    for ( int bid = 1; bid <= num_boundaryNodeIds; bid++ ) {
        char key[100];
        snprintf( key, sizeof key, "BoundaryNodeId%d", bid );
        char newKey[100];
        snprintf( newKey, sizeof newKey, "NumberOfBoundaryNodes%d", bid );
        int bndDofSize = db->getScalar<int>( newKey );
        if ( bndDofSize > 0 ) {
            auto bndDofIndices = db->getVector<int>( key );
            for ( int i = 0; i < bndDofSize; i++ ) {
                mesh->boundary_info->add_node( mesh->node_ptr( bndDofIndices[i] ), bid );
            } // end for i
        }
    } // end for bid

    for ( int bid = 1; bid <= num_boundarySideIds; bid++ ) {
        char key[100];
        snprintf( key, sizeof key, "BoundarySideId%d", bid );
        char newKey[100];
        snprintf( newKey, sizeof newKey, "NumberOfBoundarySides%d", bid );
        int bndDofSize = db->getScalar<int>( newKey );
        if ( bndDofSize > 0 ) {
            auto bndDofIndices = db->getVector<int>( key );
            for ( int i = 0; i < bndDofSize; i++ ) {
                mesh->boundary_info->add_side(
                    mesh->elem_ptr( bndDofIndices[2 * i] ), bndDofIndices[( 2 * i ) + 1], bid );
            } // end for i
        }
    } // end for bid
}


} // namespace AMP
