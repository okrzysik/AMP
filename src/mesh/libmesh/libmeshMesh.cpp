#include "AMP/mesh/libmesh/libmeshMesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/MeshElementVectorIterator.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/mesh/MultiIterator.h"
#include "AMP/mesh/libmesh/initializeLibMesh.h"
#include "AMP/mesh/libmesh/libmeshElemIterator.h"
#include "AMP/mesh/libmesh/libmeshMeshElement.h"
#include "AMP/mesh/libmesh/libmeshNodeIterator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.hpp"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include "ProfilerApp.h"


// LibMesh include
DISABLE_WARNINGS
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/exodusII_io_helper.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/parallel.h"
ENABLE_WARNINGS


namespace AMP::Mesh {


/********************************************************
 * Constructors                                          *
 ********************************************************/
libmeshMesh::libmeshMesh( std::shared_ptr<const MeshParameters> params )
    : Mesh( params ), d_pos_hash( 0 )
{
    PROFILE_START( "constructor" );
    this->d_max_gcw = 1;
    // Check for valid inputs
    AMP_INSIST( d_comm != AMP_MPI( AMP_COMM_NULL ), "Communicator must be set" );
    // Intialize libMesh
    d_libMeshComm = std::make_shared<libMesh::Parallel::Communicator>( d_comm.getCommunicator() );
    libmeshInit   = std::make_shared<initializeLibMesh>( d_comm );
    // Load the mesh
    AMP_INSIST( params.get(), "Params must not be null" );
    auto db = params->getDatabase();
    AMP_INSIST( db.get(), "Database must exist" );
    if ( db.get() ) {
        // Database exists
        AMP_INSIST( db->keyExists( "dim" ), "Variable 'dim' must be set in the database" );
        AMP_INSIST( db->keyExists( "MeshName" ), "MeshName must exist in input database" );
        PhysicalDim = db->getScalar<int>( "dim" );
        d_name      = db->getString( "MeshName" );
        AMP_INSIST( PhysicalDim > 0 && PhysicalDim < 10, "Invalid dimension" );
        GeomDim = (GeomType) PhysicalDim;
        // Create the libMesh objects
        d_libMesh = std::make_shared<libMesh::Mesh>( *d_libMeshComm, PhysicalDim );
        if ( db->keyExists( "FileName" ) ) {
            // Read an existing mesh
            d_libMesh->read( db->getString( "FileName" ) );
        } else if ( db->keyExists( "Generator" ) ) {
            // Generate a new mesh
            std::string generator = db->getString( "Generator" );
            if ( generator.compare( "cube" ) == 0 ) {
                // Generate a cube mesh
                AMP_INSIST( PhysicalDim == 3,
                            "libMesh cube generation currently supports only 3d meshes" );
                AMP_INSIST( db->keyExists( "size" ),
                            "Variable 'size' must be set in the database" );
                auto size = db->getVector<int>( "size" );
                AMP_INSIST( size.size() == (size_t) PhysicalDim,
                            "Variable 'size' must by an integer array of size dim" );
                AMP_INSIST( db->keyExists( "xmin" ),
                            "Variable 'xmin' must be set in the database" );
                auto xmin = db->getVector<double>( "xmin" );
                AMP_INSIST( xmin.size() == (size_t) PhysicalDim,
                            "Variable 'xmin' must by an integer array of size dim" );
                AMP_INSIST( db->keyExists( "xmax" ),
                            "Variable 'xmax' must be set in the database" );
                auto xmax = db->getVector<double>( "xmax" );
                AMP_INSIST( xmax.size() == (size_t) PhysicalDim,
                            "Variable 'xmax' must by an integer array of size dim" );
                libMesh::MeshTools::Generation::build_cube( *d_libMesh,
                                                            size[0],
                                                            size[1],
                                                            size[2],
                                                            xmin[0],
                                                            xmax[0],
                                                            xmin[1],
                                                            xmax[1],
                                                            xmin[2],
                                                            xmax[2],
                                                            libMesh::HEX8 );
            } else {
                AMP_ERROR( std::string( "Unknown libmesh generator: " ) + generator );
            }
        } else {
            AMP_ERROR( "Unable to construct mesh with given parameters" );
        }
        // Initialize all of the internal data
        initialize();
        // Displace the mesh
        std::vector<double> displacement( PhysicalDim, 0.0 );
        if ( db->keyExists( "x_offset" ) )
            displacement[0] = db->getScalar<double>( "x_offset" );
        if ( db->keyExists( "y_offset" ) )
            displacement[1] = db->getScalar<double>( "y_offset" );
        if ( db->keyExists( "z_offset" ) )
            displacement[2] = db->getScalar<double>( "z_offset" );
        bool test = false;
        for ( auto &elem : displacement ) {
            if ( elem != 0.0 )
                test = true;
        }
        if ( test )
            displaceMesh( displacement );
    } else {
        AMP_ERROR( "Error: params must contain a database object" );
    }
    // Get the global ranks for the comm to make sure it is set
    auto globalRanks = getComm().globalRanks();
    AMP_ASSERT( !globalRanks.empty() );
    PROFILE_STOP( "constructor" );
}
libmeshMesh::libmeshMesh( std::shared_ptr<libMesh::Mesh> mesh, const std::string &name )
    : d_pos_hash( 0 ), d_libMesh( mesh )
{
    // Set the base properties
#ifdef AMP_USE_MPI
    this->d_comm = AMP_MPI( mesh->comm().get() );
    AMP_ASSERT( !d_comm.isNull() );
#else
    this->d_comm = AMP_MPI( AMP_COMM_SELF );
#endif
    this->setMeshID();
    this->d_name      = name;
    this->d_max_gcw   = 1;
    this->PhysicalDim = d_libMesh->mesh_dimension();
    this->GeomDim     = (GeomType) PhysicalDim;
    // Initialize all of the internal data
    initialize();
}


/********************************************************
 * Destructor                                            *
 ********************************************************/
libmeshMesh::~libmeshMesh()
{
    // First we need to destroy the elements, surface sets, and boundary sets
    for ( int i = 0; i < 4; i++ ) {
        d_localElements[i].reset();
        d_ghostElements[i].reset();
    }
    d_localSurfaceElements.clear();
    d_ghostSurfaceElements.clear();
    d_boundarySets.clear();
    // We need to clear all libmesh objects before libmeshInit
    d_libMesh.reset();
    libmeshInit.reset();
}


/********************************************************
 * Return the class name                                 *
 ********************************************************/
std::string libmeshMesh::meshClass() const { return "libmeshMesh"; }


/********************************************************
 * Function to copy the mesh                             *
 ********************************************************/
std::unique_ptr<Mesh> libmeshMesh::clone() const { return std::make_unique<libmeshMesh>( *this ); }


/********************************************************
 * Function to initialize the libmeshMesh object         *
 ********************************************************/
void libmeshMesh::initialize()
{
    PROFILE_START( "initialize" );
    // Verify libmesh's rank and size agrees with the rank and size of the comm of the mesh
    AMP_INSIST( (int) d_libMesh->n_processors() == d_comm.getSize(),
                "size of the mesh does not agree with libmesh" );
    AMP_INSIST( (int) d_libMesh->processor_id() == d_comm.getRank(),
                "rank of the mesh does not agree with libmesh" );
    // Count the elements
    n_local  = std::vector<size_t>( PhysicalDim + 1, 0 );
    n_global = std::vector<size_t>( PhysicalDim + 1, 0 );
    n_ghost  = std::vector<size_t>( PhysicalDim + 1, 0 );
    for ( int i = 0; i <= (int) GeomDim; i++ ) {
        if ( i == (int) GeomType::Vertex ) {
            // We are counting the nodes
            n_local[i]  = d_libMesh->n_local_nodes();
            n_global[i] = d_libMesh->parallel_n_nodes();
            n_ghost[i] =
                std::distance( d_libMesh->nodes_begin(), d_libMesh->nodes_end() ) - n_local[i];
            auto pos = d_libMesh->nodes_begin();
            auto end = d_libMesh->nodes_end();
            int N    = 0;
            while ( pos != end ) {
                N++;
                ++pos;
            }
            n_ghost[i] = N - n_local[i];
            AMP_INSIST( n_local[i] > 0, "We currently require at least 1 node on each processor" );
        } else if ( i == (int) GeomDim ) {
            // We are counting the elements
            n_local[i]  = d_libMesh->n_local_elem();
            n_global[i] = d_libMesh->parallel_n_elem();
            n_ghost[i]  = std::distance( d_libMesh->elements_begin(), d_libMesh->elements_end() ) -
                         n_local[i];
            AMP_INSIST( n_local[i] > 0,
                        "We currently require at least 1 element on each processor" );
        } else {
            // We are counting an intermediate element (not finished)
            n_local[i]  = static_cast<size_t>( -1 );
            n_global[i] = static_cast<size_t>( -1 );
            n_ghost[i]  = static_cast<size_t>( -1 );
        }
    }
    // Compute the bounding box of the mesh
    d_box_local = std::vector<double>( PhysicalDim * 2 );
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_box_local[2 * i + 0] = 1e200;
        d_box_local[2 * i + 1] = -1e200;
    }
    auto node_pos = d_libMesh->local_nodes_begin();
    auto node_end = d_libMesh->local_nodes_end();
    while ( node_pos != node_end ) {
        libMesh::Node *node = *node_pos;
        for ( int i = 0; i < PhysicalDim; i++ ) {
            double x = ( *node )( i );
            if ( x < d_box_local[2 * i + 0] ) {
                d_box_local[2 * i + 0] = x;
            }
            if ( x > d_box_local[2 * i + 1] ) {
                d_box_local[2 * i + 1] = x;
            }
        }
        ++node_pos;
    }
    d_box = std::vector<double>( PhysicalDim * 2 );
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_box[2 * i + 0] = d_comm.minReduce( d_box_local[2 * i + 0] );
        d_box[2 * i + 1] = d_comm.maxReduce( d_box_local[2 * i + 1] );
    }
    // Construct the element neighbor information
    d_libMesh->find_neighbors();
    // Construct the node neighbor information
    neighborNodeIDs = std::vector<unsigned int>( n_local[0], (unsigned int) -1 );
    neighborNodes   = std::vector<std::vector<libMesh::Node *>>( n_local[0] );
    node_pos        = d_libMesh->local_nodes_begin();
    node_end        = d_libMesh->local_nodes_end();
    size_t i        = 0;
    while ( node_pos != node_end ) {
        auto node          = node_pos.operator*();
        neighborNodeIDs[i] = node->id();
        ++node_pos;
        i++;
    }
    AMP::Utilities::quicksort( neighborNodeIDs );
    auto elem_pos = d_libMesh->elements_begin();
    auto elem_end = d_libMesh->elements_end();
    std::vector<std::set<unsigned int>> tmpNeighborNodes( n_local[0] );
    int rank = d_comm.getRank();
    while ( elem_pos != elem_end ) {
        auto elem = elem_pos.operator*();
        for ( i = 0; i < elem->n_nodes(); i++ ) {
            auto node = elem->node_ptr( i );
            int p_id  = node->processor_id();
            if ( p_id == rank ) {
                int j = AMP::Utilities::findfirst( neighborNodeIDs, node->id() );
                AMP_ASSERT( neighborNodeIDs[j] == node->id() );
                for ( unsigned int k = 0; k < elem->n_nodes(); k++ ) {
                    if ( k == i )
                        continue;
                    libMesh::Node *node2 = elem->node_ptr( k );
                    tmpNeighborNodes[j].insert( node2->id() );
                }
            }
        }
        ++elem_pos;
    }
    for ( i = 0; i < neighborNodeIDs.size(); i++ ) {
        neighborNodes[i] = std::vector<libMesh::Node *>( tmpNeighborNodes[i].size() );
        int j            = 0;
        for ( auto it = tmpNeighborNodes[i].begin(); it != tmpNeighborNodes[i].end(); ++it ) {
            neighborNodes[i][j] = d_libMesh->node_ptr( *it );
            j++;
        }
    }
    // Construct the list of elements of type side or edge
    for ( i = 0; i <= (size_t) GeomDim; i++ ) {
        auto type = (GeomType) i;
        if ( type == GeomType::Vertex || type == GeomDim )
            continue;
        // Get a unique list of all elements of the desired type
        std::set<MeshElement> element_list;
        auto it = getIterator( GeomDim, 1 );
        for ( size_t j = 0; j < it.size(); j++ ) {
            auto tmp = it->getElements( type );
            for ( auto &elem : tmp )
                element_list.insert( elem );
            ++it;
        }
        // Split the new elements into the local and ghost lists
        size_t N_local = 0;
        size_t N_ghost = 0;
        for ( auto elem : element_list ) {
            MeshElementID id = elem.globalID();
            if ( id.is_local() )
                N_local++;
            else
                N_ghost++;
        }
        size_t N_global = d_comm.sumReduce( N_local );
        AMP_ASSERT( N_global >= n_global[static_cast<int>( GeomDim )] );
        auto local_elements = std::make_shared<std::vector<MeshElement>>( N_local );
        auto ghost_elements = std::make_shared<std::vector<MeshElement>>( N_ghost );
        N_local             = 0;
        N_ghost             = 0;
        for ( const auto &elem : element_list ) {
            MeshElementID id = elem.globalID();
            if ( id.is_local() ) {
                local_elements->operator[]( N_local ) = elem;
                N_local++;
            } else {
                ghost_elements->operator[]( N_ghost ) = elem;
                N_ghost++;
            }
        }
        AMP::Utilities::quicksort( *local_elements ); // Sort elements for searching
        AMP::Utilities::quicksort( *ghost_elements ); // Sort elements for searching
        d_localElements[i] = local_elements;
        d_ghostElements[i] = ghost_elements;
        n_local[i]         = local_elements->size();
        n_global[i]        = d_comm.sumReduce( n_local[i] );
        n_ghost[i]         = ghost_elements->size();
    }
    /*for (int i=0; i<=(int)GeomDim; i++) {
        for (int j=0; j<d_comm.getSize(); j++) {
            d_comm.barrier();
            if ( d_comm.getRank()==j )
                printf("%i, %i, %i, %i\n",i,n_local[i],n_ghost[i],n_global[i]);
            d_comm.barrier();
            std::cout << "";
        }
    }*/
    // Construct the boundary elements for Node and Elem
    d_localSurfaceElements =
        std::vector<std::shared_ptr<std::vector<MeshElement>>>( (int) GeomDim + 1 );
    d_ghostSurfaceElements =
        std::vector<std::shared_ptr<std::vector<MeshElement>>>( (int) GeomDim + 1 );
    elem_pos = d_libMesh->elements_begin();
    elem_end = d_libMesh->elements_end();
    std::set<libMesh::Elem *> localBoundaryElements;
    std::set<libMesh::Elem *> ghostBoundaryElements;
    std::set<libMesh::Node *> localBoundaryNodes;
    std::set<libMesh::Node *> ghostBoundaryNodes;
    while ( elem_pos != elem_end ) {
        auto element = *elem_pos;
        if ( element->on_boundary() ) {
            if ( (int) element->processor_id() == rank )
                localBoundaryElements.insert( element );
            else
                ghostBoundaryElements.insert( element );
            for ( unsigned int si = 0; si < element->n_sides(); si++ ) {
                if ( element->neighbor_ptr( si ) == nullptr ) {
                    auto side = element->build_side_ptr( si );
                    for ( unsigned int j = 0; j < side->n_nodes(); j++ ) {
                        auto node = side->node_ptr( j );
                        if ( (int) node->processor_id() == rank )
                            localBoundaryNodes.insert( node );
                        else
                            ghostBoundaryNodes.insert( node );
                    }
                }
            }
        }
        elem_pos++;
    }
    auto GeomDim2 = static_cast<int>( GeomDim );
    d_localSurfaceElements[GeomDim2] =
        std::make_shared<std::vector<MeshElement>>( localBoundaryElements.size() );
    auto elem_iterator = localBoundaryElements.begin();
    for ( i = 0; i < localBoundaryElements.size(); i++ ) {
        ( *d_localSurfaceElements[GeomDim2] )[i] = libmeshMeshElement(
            PhysicalDim, GeomDim, (void *) *elem_iterator, rank, d_meshID, this );
        ++elem_iterator;
    }
    AMP::Utilities::quicksort( *d_localSurfaceElements[GeomDim2] );
    d_ghostSurfaceElements[GeomDim2] =
        std::make_shared<std::vector<MeshElement>>( ghostBoundaryElements.size() );
    elem_iterator = ghostBoundaryElements.begin();
    for ( i = 0; i < ghostBoundaryElements.size(); i++ ) {
        ( *d_ghostSurfaceElements[GeomDim2] )[i] = libmeshMeshElement(
            PhysicalDim, GeomDim, (void *) *elem_iterator, rank, d_meshID, this );
        ++elem_iterator;
    }
    AMP::Utilities::quicksort( *d_ghostSurfaceElements[GeomDim2] );
    d_localSurfaceElements[0] =
        std::make_shared<std::vector<MeshElement>>( localBoundaryNodes.size() );
    auto node_iterator = localBoundaryNodes.begin();
    for ( i = 0; i < localBoundaryNodes.size(); i++ ) {
        ( *d_localSurfaceElements[0] )[i] = libmeshMeshElement(
            PhysicalDim, GeomType::Vertex, (void *) *node_iterator, rank, d_meshID, this );
        ++node_iterator;
    }
    AMP::Utilities::quicksort( *d_localSurfaceElements[0] );
    d_ghostSurfaceElements[0] =
        std::make_shared<std::vector<MeshElement>>( ghostBoundaryNodes.size() );
    node_iterator = ghostBoundaryNodes.begin();
    for ( i = 0; i < ghostBoundaryNodes.size(); i++ ) {
        ( *d_ghostSurfaceElements[0] )[i] = libmeshMeshElement(
            PhysicalDim, GeomType::Vertex, (void *) *node_iterator, rank, d_meshID, this );
        ++node_iterator;
    }
    AMP::Utilities::quicksort( *d_ghostSurfaceElements[0] );
    // Construct the boundary elements for all other types
    // An face or edge is on the boundary if all of its nodes are on the surface
    size_t element_surface_global_size =
        d_comm.sumReduce( d_localSurfaceElements[GeomDim2]->size() );
    for ( int type2 = 1; type2 < (int) GeomDim; type2++ ) {
        auto type = (GeomType) type2;
        std::set<MeshElement> local, ghost;
        MeshIterator it = getIterator( type, 0 );
        for ( i = 0; i < it.size(); i++ ) {
            auto nodes = it->getElements( GeomType::Vertex );
            AMP_ASSERT( !nodes.empty() );
            bool on_boundary = true;
            for ( auto &node : nodes ) {
                if ( !node.isOnSurface() )
                    on_boundary = false;
            }
            if ( on_boundary ) {
                if ( it->globalID().is_local() )
                    local.insert( *it );
                else
                    ghost.insert( *it );
            }
            ++it;
        }
        d_localSurfaceElements[type2] =
            std::make_shared<std::vector<MeshElement>>( local.begin(), local.end() );
        d_ghostSurfaceElements[type2] =
            std::make_shared<std::vector<MeshElement>>( ghost.begin(), ghost.end() );
        AMP::Utilities::quicksort( *d_localSurfaceElements[type2] );
        AMP::Utilities::quicksort( *d_ghostSurfaceElements[type2] );
        size_t local_size  = d_localSurfaceElements[type2]->size();
        size_t global_size = d_comm.sumReduce( local_size );
        AMP_ASSERT( global_size >= element_surface_global_size );
    }
    // Construct the boundary lists
    auto libmesh_bids = d_libMesh->boundary_info->get_boundary_ids();
    std::vector<short int> bids( libmesh_bids.begin(), libmesh_bids.end() );
    Utilities::quicksort( bids );
    std::vector<short int> side_ids;
    std::vector<short int> node_ids;
    d_libMesh->boundary_info->build_side_boundary_ids( side_ids );
    d_libMesh->boundary_info->build_node_list_from_side_list();
    d_libMesh->boundary_info->build_node_boundary_ids( node_ids );
    for ( int type2 = 0; type2 <= (int) GeomDim; type2++ ) {
        auto type     = (GeomType) type2;
        auto iterator = getIterator( type, 0 );
        for ( auto &bid : bids ) {
            auto id = (int) bid;
            // Count the number of elements on the given boundary
            auto curElem = iterator.begin();
            auto endElem = iterator.end();
            int N        = 0;
            while ( curElem != endElem ) {
                if ( curElem->isOnBoundary( id ) )
                    N++;
                ++curElem;
            }
            // Create the boundary list
            auto list = std::make_shared<std::vector<MeshElement>>( N );
            curElem   = iterator.begin();
            endElem   = iterator.end();
            N         = 0;
            while ( curElem != endElem ) {
                if ( curElem->isOnBoundary( id ) ) {
                    list->operator[]( N ) = *curElem;
                    N++;
                }
                ++curElem;
            }
            // Store the list
            auto mapid = std::pair<int, GeomType>( id, type );
            auto entry = std::make_pair( mapid, list );
            d_boundarySets.insert( entry );
        }
    }
    // Get a list of all block ids
    std::set<int> block_ids;
    elem_pos = d_libMesh->elements_begin();
    elem_end = d_libMesh->elements_end();
    while ( elem_pos != elem_end ) {
        libMesh::Elem *element = *elem_pos;
        int id                 = element->subdomain_id();
        block_ids.insert( id );
        ++elem_pos;
    }
    std::vector<int> send_list( block_ids.begin(), block_ids.end() );
    size_t recv_size = d_comm.sumReduce( send_list.size() );
    std::vector<int> recv_list( recv_size, 0 );
    d_comm.allGather( &send_list[0], send_list.size(), &recv_list[0] );
    for ( auto &elem : recv_list )
        block_ids.insert( elem );
    d_block_ids = std::vector<int>( block_ids.begin(), block_ids.end() );
    PROFILE_STOP( "initialize" );
    // Test getting some common iterators
    auto it1 = getIterator( AMP::Mesh::GeomType::Vertex, 0 );
    auto it2 = getIterator( this->GeomDim, 0 );
    AMP_ASSERT( it1.size() == n_local[0] );
    AMP_ASSERT( it2.size() == n_local[static_cast<int>( this->GeomDim )] );
}


/********************************************************
 * Function to estimate the mesh size                    *
 ********************************************************/
size_t libmeshMesh::estimateMeshSize( std::shared_ptr<const MeshParameters> params )
{
    auto database = params->getDatabase();
    AMP_ASSERT( database );
    size_t NumberOfElements = 0;
    if ( database->keyExists( "NumberOfElements" ) ) {
        // User specified the number of elements, this should override everything
        NumberOfElements = (size_t) database->getScalar<int>( "NumberOfElements" );
    } else if ( database->keyExists( "FileName" ) ) {
        // Read an existing mesh
        auto fname = database->getString( "FileName" );
        if ( fname.rfind( ".exd" ) < fname.size() || fname.rfind( ".e" ) < fname.size() ) {
            libMesh::Parallel::Communicator comm;
            libMesh::ExodusII_IO_Helper exio_helper( comm, false, true );
            exio_helper.open( fname.c_str(), true ); // Open the exodus file, if possible
#if ( LIBMESH_MINOR_VERSION < 6 )
            exio_helper.read_header(); // Read the header
#else
            // the read_header function above now does something different!!
            exio_helper.read_and_store_header_info(); // Read the header
#endif
            exio_helper.close(); // Close the file
            NumberOfElements = exio_helper.num_elem;
            AMP_ASSERT( NumberOfElements > 0 );
        } else {
            AMP_ERROR( "Unkown mesh type, use key NumberOfElements to specify the mesh size" );
        }
    } else if ( database->keyExists( "Generator" ) ) {
        // Generate a new mesh
        auto generator = database->getString( "Generator" );
        if ( generator.compare( "cube" ) == 0 ) {
            // Generate a cube mesh
            AMP_INSIST( database->keyExists( "size" ),
                        "Variable 'size' must be set in the database" );
            auto size        = database->getVector<int>( "size" );
            NumberOfElements = 1;
            for ( auto &elem : size )
                NumberOfElements *= elem;
        } else {
            AMP_ERROR( std::string( "Unknown libmesh generator: " ) + generator );
        }
    } else {
        AMP_ERROR( "Unable to construct mesh with given parameters" );
    }
    // Adjust the number of elements by a weight if desired
    if ( database->keyExists( "Weight" ) ) {
        auto weight      = database->getScalar<double>( "Weight" );
        NumberOfElements = (size_t) ceil( weight * ( (double) NumberOfElements ) );
    }
    return NumberOfElements;
}


/****************************************************************
 * Estimate the maximum number of processors                     *
 ****************************************************************/
size_t libmeshMesh::maxProcs( std::shared_ptr<const MeshParameters> params )
{
    return estimateMeshSize( params );
}


/********************************************************
 * Return the number of elements                         *
 ********************************************************/
size_t libmeshMesh::numLocalElements( const GeomType type ) const
{
    auto n = n_local[static_cast<int>( type )];
    if ( n == static_cast<size_t>( -1 ) )
        AMP_ERROR( "numLocalElements is not implemented for this type" );
    return n;
}
size_t libmeshMesh::numGlobalElements( const GeomType type ) const
{
    auto n = n_global[static_cast<int>( type )];
    if ( n == static_cast<size_t>( -1 ) )
        AMP_ERROR( "numLocalElements is not implemented for this type" );
    return n;
}
size_t libmeshMesh::numGhostElements( const GeomType type, int gcw ) const
{
    if ( gcw == 0 )
        return 0;
    if ( gcw > 1 )
        AMP_ERROR( "Libmesh only supports a ghost cell width of 1" );
    auto n = n_ghost[static_cast<int>( type )];
    if ( n == static_cast<size_t>( -1 ) )
        AMP_ERROR( "numLocalElements is not implemented for this type" );
    return n;
}


/********************************************************
 * Return an iterator over the given geometric type      *
 ********************************************************/
MeshIterator libmeshMesh::getIterator( const GeomType type, const int gcw ) const
{
    MeshIterator it;
    int i = static_cast<int>( type );
    if ( static_cast<int>( type ) == PhysicalDim ) {
        // This is a libMesh element
        if ( gcw == 0 ) {
            auto begin = d_libMesh->local_elements_begin();
            auto end   = d_libMesh->local_elements_end();
            it         = libmeshElemIterator( this, gcw, begin, end, begin, n_local[i], 0 );
        } else if ( gcw == 1 ) {
            auto begin = d_libMesh->elements_begin();
            auto end   = d_libMesh->elements_end();
            it = libmeshElemIterator( this, gcw, begin, end, begin, n_local[i] + n_ghost[i], 0 );
        } else {
            AMP_ERROR( "Unsupported ghost cell width" );
        }
    } else if ( type == GeomType::Vertex ) {
        // This is a libMesh node
        if ( gcw == 0 ) {
            auto begin = d_libMesh->local_nodes_begin();
            auto end   = d_libMesh->local_nodes_end();
            it         = libmeshNodeIterator( this, gcw, begin, end, begin, n_local[i], 0 );
        } else if ( gcw == 1 ) {
            auto begin = d_libMesh->nodes_begin();
            auto end   = d_libMesh->nodes_end();
            it = libmeshNodeIterator( this, gcw, begin, end, begin, n_local[i] + n_ghost[i], 0 );
        } else {
            AMP_ERROR( "Unsupported ghost cell width" );
        }
    } else {
        // All other types require a pre-constructed list
        if ( gcw == 0 ) {
            it = MultiVectorIterator( d_localElements[i], 0 );
        } else if ( gcw == 1 ) {
            std::vector<MeshIterator> iterators( 2 );
            iterators[0] = MultiVectorIterator( d_localElements[i], 0 );
            iterators[1] = MultiVectorIterator( d_ghostElements[i], 0 );
            it           = MultiIterator( iterators, 0 );
        } else {
            AMP_ERROR( "Unsupported ghost cell width" );
        }
    }
    return it;
}


/********************************************************
 * Return an iterator over the given boundary ids        *
 * Note: we have not programmed this for ghosts yet      *
 ********************************************************/
MeshIterator libmeshMesh::getSurfaceIterator( const GeomType type, const int gcw ) const
{
    AMP_ASSERT( type <= GeomDim );
    auto local = d_localSurfaceElements[(int) type];
    auto ghost = d_ghostSurfaceElements[(int) type];
    if ( local.get() == nullptr || ghost.get() == nullptr )
        AMP_ERROR( "Surface iterator over the given geometry type is not supported" );
    if ( gcw == 0 ) {
        return MultiVectorIterator( local, 0 );
    } else if ( gcw == 1 ) {
        std::vector<MeshIterator> iterators( 2 );
        iterators[0] = MultiVectorIterator( local, 0 );
        iterators[1] = MultiVectorIterator( ghost, 0 );
        return MultiIterator( iterators, 0 );
    } else {
        AMP_ERROR( "libmesh has maximum ghost width of 1" );
    }
    return MeshIterator();
}


/********************************************************
 * Return an iterator over the given boundary ids        *
 * Note: we have not programmed this for ghosts yet      *
 ********************************************************/
std::vector<int> libmeshMesh::getBoundaryIDs() const
{
    auto libmesh_bids = d_libMesh->boundary_info->get_boundary_ids();
    std::vector<int> bids( libmesh_bids.size(), 0 );
    auto it = libmesh_bids.begin();
    for ( auto &bid : bids ) {
        bid = *it;
        ++it;
    }
    return bids;
}
MeshIterator
libmeshMesh::getBoundaryIDIterator( const GeomType type, const int id, const int gcw ) const
{
    AMP_INSIST( gcw == 0, "Iterator over ghost boundary elements is not supported yet" );
    auto mapid = std::pair<int, GeomType>( id, type );
    auto list  = std::make_shared<std::vector<MeshElement>>();
    auto it    = d_boundarySets.find( mapid );
    if ( it != d_boundarySets.end() )
        list = it->second;
    return MultiVectorIterator( list, 0 );
}


/********************************************************
 * Return an iterator over the given block ids           *
 ********************************************************/
std::vector<int> libmeshMesh::getBlockIDs() const { return d_block_ids; }
MeshIterator
libmeshMesh::getBlockIDIterator( const GeomType type, const int id, const int gcw ) const
{
    if ( d_block_ids.size() == 1 ) {
        if ( d_block_ids[0] == id )
            return getIterator( type, gcw );
        else
            return MeshIterator();
    } else {
        AMP_ERROR( "getBlockIDIterator is not implemented yet" );
    }
    return MeshIterator();
}


/********************************************************
 * Return pointers to the neighbor nodes give a node id  *
 ********************************************************/
std::vector<libMesh::Node *> libmeshMesh::getNeighborNodes( const MeshElementID &id ) const
{
    AMP_INSIST( id.type() == GeomType::Vertex, "This function is for nodes" );
    AMP_INSIST( id.meshID() == d_meshID, "Unknown mesh" );
    AMP_INSIST( id.is_local(), "Only owned nodes can return their neighbor lists" );
    int i = AMP::Utilities::findfirst( neighborNodeIDs, id.local_id() );
    AMP_ASSERT( neighborNodeIDs[i] == id.local_id() );
    return neighborNodes[i];
}


/********************************************************
 * Function to return the element given an ID            *
 ********************************************************/
MeshElement libmeshMesh::getElement( const MeshElementID &elem_id ) const
{
    MeshID mesh_id = elem_id.meshID();
    AMP_INSIST( mesh_id == d_meshID, "mesh id must match the mesh id of the element" );
    unsigned int rank = d_comm.getRank();
    if ( (int) elem_id.type() == PhysicalDim ) {
        // This is a libMesh element
        auto element = d_libMesh->elem_ptr( elem_id.local_id() );
        return libmeshMeshElement(
            PhysicalDim, elem_id.type(), (void *) element, rank, mesh_id, this );
    } else if ( elem_id.type() == GeomType::Vertex ) {
        // This is a libMesh node
        auto node = d_libMesh->node_ptr( elem_id.local_id() );
        return libmeshMeshElement(
            PhysicalDim, elem_id.type(), (void *) node, rank, mesh_id, this );
    }
    // All other types are stored in sorted lists
    std::shared_ptr<std::vector<MeshElement>> list;
    if ( (int) elem_id.owner_rank() == d_comm.getRank() )
        list = d_localElements[static_cast<int>( elem_id.type() )];
    else
        list = d_ghostElements[static_cast<int>( elem_id.type() )];
    size_t n = list->size();
    AMP_ASSERT( n > 0 );
    auto x = &( list->operator[]( 0 ) ); // Use the pointer for speed
    if ( x[0] == elem_id )
        return x[0];
    size_t lower = 0;
    size_t upper = n - 1;
    size_t index;
    while ( ( upper - lower ) != 1 ) {
        index = ( upper + lower ) / 2;
        if ( x[index] >= elem_id )
            upper = index;
        else
            lower = index;
    }
    index = upper;
    if ( x[index] == elem_id )
        return x[index];
    if ( elem_id.is_local() )
        AMP_ERROR( "Local element not found" );
    return MeshElement();
}


/********************************************************
 * Displace a mesh                                       *
 ********************************************************/
Mesh::Movable libmeshMesh::isMeshMovable() const { return Mesh::Movable::Displace; }
uint64_t libmeshMesh::positionHash() const { return d_pos_hash; }
void libmeshMesh::displaceMesh( const std::vector<double> &x_in )
{
    // Check x
    AMP_INSIST( (short int) x_in.size() == PhysicalDim,
                "Displacement vector size should match PhysicalDim" );
    auto x = x_in;
    d_comm.minReduce( &x[0], x.size() );
    for ( size_t i = 0; i < x.size(); i++ )
        AMP_INSIST( fabs( x[i] - x_in[i] ) < 1e-12, "x does not match on all processors" );
    // Move the mesh
    auto cur = d_libMesh->nodes_begin();
    auto end = d_libMesh->nodes_end();
    while ( cur != end ) {
        auto d_Node = *cur;
        for ( size_t i = 0; i < x.size(); i++ )
            ( *d_Node )( i ) += x[i];
        ++cur;
    }
    // Update the bounding box
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_box[2 * i + 0] += x[i];
        d_box[2 * i + 1] += x[i];
        d_box_local[2 * i + 0] += x[i];
        d_box_local[2 * i + 1] += x[i];
    }
    d_pos_hash++;
}
void libmeshMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr x )
{
    // Create the position vector with the necessary ghost nodes
    auto DOFs = AMP::Discretization::simpleDOFManager::create(
        shared_from_this(),
        getIterator( AMP::Mesh::GeomType::Vertex, 1 ),
        getIterator( AMP::Mesh::GeomType::Vertex, 0 ),
        PhysicalDim );
    auto nodalVariable = std::make_shared<AMP::LinearAlgebra::Variable>( "tmp_pos" );
    auto displacement  = AMP::LinearAlgebra::createVector( DOFs, nodalVariable, false );
    std::vector<size_t> dofs1( PhysicalDim );
    std::vector<size_t> dofs2( PhysicalDim );
    auto cur  = getIterator( AMP::Mesh::GeomType::Vertex, 0 );
    auto end  = cur.end();
    auto DOFx = x->getDOFManager();
    std::vector<double> data( PhysicalDim );
    while ( cur != end ) {
        auto id = cur->globalID();
        DOFx->getDOFs( id, dofs1 );
        DOFs->getDOFs( id, dofs2 );
        x->getValuesByGlobalID( PhysicalDim, &dofs1[0], &data[0] );
        displacement->setValuesByGlobalID( PhysicalDim, &dofs2[0], &data[0] );
        ++cur;
    }
    displacement->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    // Move all nodes (including the ghost nodes)
    auto node_cur = d_libMesh->nodes_begin();
    auto node_end = d_libMesh->nodes_end();
    int rank      = d_comm.getRank();
    while ( node_cur != node_end ) {
        auto node = *node_cur;
        // Create the element id
        auto owner_rank = node->processor_id();
        auto local_id   = node->id();
        bool is_local   = (int) owner_rank == rank;
        AMP::Mesh::MeshElementID id(
            is_local, AMP::Mesh::GeomType::Vertex, local_id, owner_rank, d_meshID );
        // Get the position of the point
        DOFs->getDOFs( id, dofs2 );
        displacement->getValuesByGlobalID( PhysicalDim, &dofs2[0], &data[0] );
        // Move the point
        for ( int i = 0; i < PhysicalDim; i++ )
            ( *node )( i ) += data[i];
        ++node_cur;
    }
    // Compute the bounding box of the mesh
    d_box_local = std::vector<double>( PhysicalDim * 2 );
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_box_local[2 * i + 0] = 1e200;
        d_box_local[2 * i + 1] = -1e200;
    }
    node_cur = d_libMesh->local_nodes_begin();
    node_end = d_libMesh->local_nodes_end();
    while ( node_cur != node_end ) {
        libMesh::Node *node = *node_cur;
        for ( int i = 0; i < PhysicalDim; i++ ) {
            double x = ( *node )( i );
            if ( x < d_box_local[2 * i + 0] ) {
                d_box_local[2 * i + 0] = x;
            }
            if ( x > d_box_local[2 * i + 1] ) {
                d_box_local[2 * i + 1] = x;
            }
        }
        ++node_cur;
    }
    d_box = std::vector<double>( PhysicalDim * 2 );
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_box[2 * i + 0] = d_comm.minReduce( d_box_local[2 * i + 0] );
        d_box[2 * i + 1] = d_comm.maxReduce( d_box_local[2 * i + 1] );
    }
    d_pos_hash++;
}


/****************************************************************
 * Check if two meshes are equal                                 *
 ****************************************************************/
bool libmeshMesh::operator==( const Mesh &rhs ) const
{
    // Check if &rhs == this
    if ( this == &rhs )
        return true;
    // Check if we can cast to a MultiMesh
    auto mesh = dynamic_cast<const libmeshMesh *>( &rhs );
    if ( !mesh )
        return false;
    // Perform comparison on sub-meshes
    AMP_ERROR( "Not finished" );
    return false;
}


/****************************************************************
 * Write restart data                                            *
 ****************************************************************/
void libmeshMesh::writeRestart( int64_t ) const
{
    AMP_ERROR( "writeRestart is not implimented for libmeshMesh" );
}


} // namespace AMP::Mesh
