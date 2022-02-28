#include "ReadTestMesh.h"

// LibMesh include
DISABLE_WARNINGS
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
    FILE *fp = fopen( mesh_file.c_str(), "r" );
    char str[256];
    int n;

    // Mesh {
    n = fscanf( fp, "%s {", str );
    AMP_INSIST( ( n == 1 ), "Error while reading the file" );

    int num_nodes;
    // NumberOfNodes = <val>
    n = fscanf( fp, "%s = %d", str, &num_nodes );

    mesh->reserve_nodes( num_nodes );

    for ( int i = 0; i < num_nodes; i++ ) {
        double point[3];
        // PointK = x, y, z
        n = fscanf( fp, "%s = %lf, %lf, %lf", str, &( point[0] ), &( point[1] ), &( point[2] ) );
        mesh->add_point( libMesh::Point( point[0], point[1], point[2] ), i );
    } // end for i

    int num_elem;
    // NumberOfElements = <val>
    n = fscanf( fp, "%s = %d", str, &num_elem );

    mesh->reserve_elem( num_elem );

    std::vector<std::vector<int>> elemNodeMap;

    for ( int i = 0; i < num_elem; i++ ) {
        std::vector<int> nodesForElem( 8 );
        // ElemK = i0, i1, i2, i3, i4, i5, i6, i7
        n = fscanf( fp, "%s = %d,", str, &( nodesForElem[0] ) );
        for ( int j = 1; j < 7; j++ ) {
            n = fscanf( fp, "%d,", &( nodesForElem[j] ) );
            AMP_INSIST( ( n == 1 ), "Error while reading the file" );
        } // end for j
        n = fscanf( fp, "%d", &( nodesForElem[7] ) );
        AMP_INSIST( ( n == 1 ), "Error while reading the file" );
        elemNodeMap.push_back( nodesForElem );
    } // end for i

    for ( int i = 0; i < num_elem; i++ ) {
        libMesh::Elem *newElem = new libMesh::Hex8;
        newElem->set_id( i );
        libMesh::Elem *elem = mesh->add_elem( newElem );
        for ( int j = 0; j < 8; j++ ) {
            elem->set_node( j ) = mesh->node_ptr( elemNodeMap[i][j] );
        } // end for j
    }     // end for i

    int num_boundaryNodeIds;
    // NumberOfBoundaryNodeIds = <val>
    n = fscanf( fp, "%s = %d", str, &num_boundaryNodeIds );

    for ( int bid = 1; bid <= num_boundaryNodeIds; bid++ ) {
        int bndDofSize;
        // NumberOfBoundaryNodesK = <val>
        n = fscanf( fp, "%s = %d", str, &bndDofSize );
        if ( bndDofSize > 0 ) {
            // BoundaryNodeIdK =
            n = fscanf( fp, "%s =", str );
            AMP_INSIST( ( n == 1 ), "Error while reading the file" );
            int idx;
            for ( int i = 0; i < ( bndDofSize - 1 ); i++ ) {
                n = fscanf( fp, "%d,", &idx );
                AMP_INSIST( ( n == 1 ), "Error while reading the file" );
                mesh->boundary_info->add_node( mesh->node_ptr( idx ), bid );
            } // end for i
            n = fscanf( fp, "%d", &idx );
            mesh->boundary_info->add_node( mesh->node_ptr( idx ), bid );
        }
    } // end for bid

    int num_boundarySideIds;
    // NumberOfBoundarySideIds = <val>
    n = fscanf( fp, "%s = %d", str, &num_boundarySideIds );

    for ( int bid = 1; bid <= num_boundarySideIds; bid++ ) {
        int bndDofSize;
        // NumberOfBoundarySidesK = <val>
        n = fscanf( fp, "%s = %d", str, &bndDofSize );
        if ( bndDofSize > 0 ) {
            // BoundarySideIdK =
            n = fscanf( fp, "%s =", str );
            int idxE;
            int idxS;
            for ( int i = 0; i < ( bndDofSize - 1 ); i++ ) {
                n = fscanf( fp, "%d,", &idxE );
                AMP_INSIST( ( n == 1 ), "Error while reading the file" );
                n = fscanf( fp, "%d,", &idxS );
                AMP_INSIST( ( n == 1 ), "Error while reading the file" );
                mesh->boundary_info->add_side( mesh->elem_ptr( idxE ), idxS, bid );
            } // end for i
            n = fscanf( fp, "%d,", &idxE );
            AMP_INSIST( ( n == 1 ), "Error while reading the file" );
            n = fscanf( fp, "%d", &idxS );
            AMP_INSIST( ( n == 1 ), "Error while reading the file" );
            mesh->boundary_info->add_side( mesh->elem_ptr( idxE ), idxS, bid );
        }
    } // end for bid

    fclose( fp );
}

void readTestMesh( std::shared_ptr<AMP::Database> mesh_file_db,
                   std::shared_ptr<libMesh::Mesh> mesh )
{
    auto mesh_db            = mesh_file_db->getDatabase( "Mesh" );
    int num_elem            = mesh_db->getScalar<int>( "NumberOfElements" );
    int num_nodes           = mesh_db->getScalar<int>( "NumberOfNodes" );
    int num_boundaryNodeIds = mesh_db->getScalar<int>( "NumberOfBoundaryNodeIds" );
    int num_boundarySideIds = mesh_db->getScalar<int>( "NumberOfBoundarySideIds" );

    mesh->reserve_elem( num_elem );
    mesh->reserve_nodes( num_nodes );

    for ( int i = 0; i < num_nodes; i++ ) {
        char key[100];
        snprintf( key, sizeof key, "Point%d", i );
        auto point = mesh_db->getVector<double>( key );
        mesh->add_point( libMesh::Point( point[0], point[1], point[2] ), i );
    } // end for i

    std::vector<std::vector<int>> elemNodeMap;

    for ( int i = 0; i < num_elem; i++ ) {
        char key[100];
        snprintf( key, sizeof key, "Elem%d", i );
        elemNodeMap.push_back( mesh_db->getVector<int>( key ) );
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
        int bndDofSize = mesh_db->getScalar<int>( newKey );
        if ( bndDofSize > 0 ) {
            auto bndDofIndices = mesh_db->getVector<int>( key );
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
        int bndDofSize = mesh_db->getScalar<int>( newKey );
        if ( bndDofSize > 0 ) {
            auto bndDofIndices = mesh_db->getVector<int>( key );
            for ( int i = 0; i < bndDofSize; i++ ) {
                mesh->boundary_info->add_side(
                    mesh->elem_ptr( bndDofIndices[2 * i] ), bndDofIndices[( 2 * i ) + 1], bid );
            } // end for i
        }
    } // end for bid
}
} // namespace AMP
