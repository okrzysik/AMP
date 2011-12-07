#include "ampmesh/libmesh/initializeLibMesh.h"
#include "ampmesh/libmesh/libMesh.h"
#include "ampmesh/libmesh/libMeshIterator.h"
#include "ampmesh/MeshElementVectorIterator.h"
#include "utils/MemoryDatabase.h"
#include "utils/AMPManager.h"
#include "utils/Utilities.h"


// LibMesh include
#include "mesh.h"
#include "mesh_data.h"
#include "mesh_generation.h"
#include "boundary_info.h"


namespace AMP {
namespace Mesh {


/********************************************************
* Constructors                                          *
********************************************************/
libMesh::libMesh( const MeshParameters::shared_ptr &params_in ):
    Mesh(params_in)
{
    // Check for valid inputs
    AMP_INSIST(params.get(),"Params must not be null");
    AMP_INSIST(d_comm!=AMP_MPI(AMP_COMM_NULL),"Communicator must be set");
    // Intialize libMesh (this needs to be moved out of AMPManager)
    libmeshInit = boost::shared_ptr<initializeLibMesh>(new initializeLibMesh(d_comm));
    // Load the mesh
    if ( d_db.get() ) {
        // Database exists
        AMP_INSIST(d_db->keyExists("dim"),"Variable 'dim' must be set in the database");
        AMP_INSIST(d_db->keyExists("MeshName"),"MeshName must exist in input database");
        PhysicalDim = d_db->getInteger("dim");
        d_name = d_db->getString("MeshName");
        AMP_INSIST(PhysicalDim>0&&PhysicalDim<10,"Invalid dimension");
        GeomDim = (GeomType) PhysicalDim;
        // Create the libMesh objects
        d_libMesh = boost::shared_ptr< ::Mesh>( new ::Mesh(PhysicalDim) );
        d_libMeshData = boost::shared_ptr< ::MeshData>( new ::MeshData(*d_libMesh) );
        if ( d_db->keyExists("FileName") ) {
            // Read an existing mesh
            d_libMesh->read(d_db->getString("FileName"));
        } else if ( d_db->keyExists("Generator") ) {
            // Generate a new mesh
            std::string generator = d_db->getString("Generator");
            if ( generator.compare("cube")==0 ) {
                // Generate a cube mesh
                AMP_INSIST(PhysicalDim==3,"libMesh cube generation currently supports only 3d meshes");
                AMP_INSIST(d_db->keyExists("size"),"Variable 'size' must be set in the database");
                std::vector<int> size = d_db->getIntegerArray("size");
                AMP_INSIST(size.size()==(size_t)PhysicalDim,"Variable 'size' must by an integer array of size dim");
                MeshTools::Generation::build_cube( *d_libMesh, size[0], size[1], size[2], -1., 1, -1, 1, -1, 1, HEX8 );
            } else {
                AMP_ERROR(std::string("Unknown libmesh generator: ")+generator);
            }
        } else {
            AMP_ERROR("Unable to construct mesh with given parameters");
        }
        // Initialize all of the internal data
        initialize();
        // Displace the mesh
        std::vector<double> displacement(PhysicalDim,0.0);
        if ( d_db->keyExists("x_offset") )
            displacement[0] = d_db->getDouble("x_offset");
        if ( d_db->keyExists("y_offset") )
            displacement[1] = d_db->getDouble("y_offset");
        if ( d_db->keyExists("z_offset") )
            displacement[2] = d_db->getDouble("z_offset");
        bool test = false;
        for (size_t i=0; i<displacement.size(); i++) {
            if ( displacement[i] != 0.0 )
                test = true;
        }        
        if ( test )
            displaceMesh(displacement);
    } else {
        AMP_ERROR("Error: params must contain a database object");
    }
}


/********************************************************
* De-constructor                                        *
********************************************************/
libMesh::~libMesh()
{
    // First we need to destroy the boundary sets
    d_boundarySets.clear();
    // We need to clear all libmesh objects before libmeshInit
    d_libMeshData.reset();
    d_libMesh.reset();
    libmeshInit.reset();
}


/********************************************************
* Function to copy the mesh                             *
********************************************************/
Mesh libMesh::copy() const
{
    return libMesh(*this);
}


/********************************************************
* Function to initialize the libMesh object             *
********************************************************/
void libMesh::initialize()
{
    // Verify libmesh's rank and size agrees with the rank and size of the comm of the mesh
    AMP_INSIST((int)d_libMesh->processor_id()==d_comm.getRank(),"rank of the mesh does not agree with libmesh");
    AMP_INSIST((int)d_libMesh->n_processors()==d_comm.getSize(),"size of the mesh does not agree with libmesh");
    // Count the elements 
    n_local = std::vector<size_t>(PhysicalDim+1,0);
    n_global = std::vector<size_t>(PhysicalDim+1,0);
    n_ghost = std::vector<size_t>(PhysicalDim+1,0);
    for (int i=0; i<=(int)GeomDim; i++) {
        if ( i == (int) Vertex ) {
            // We are counting the nodes
            n_local[i] = d_libMesh->n_local_nodes();
            n_global[i] = d_libMesh->parallel_n_nodes();
            ::Mesh::node_iterator pos = d_libMesh->nodes_begin();
            ::Mesh::node_iterator end = d_libMesh->nodes_end();
            int N = 0;
            while ( pos != end ) {
                N++;
                ++pos;
            }
            n_ghost[i] = N - n_local[i];
        } else if ( i == (int) GeomDim ) {
            // We are counting the elements
            n_local[i] = d_libMesh->n_local_elem();
            n_global[i] = d_libMesh->parallel_n_elem();
            ::Mesh::element_iterator pos = d_libMesh->elements_begin();
            ::Mesh::element_iterator end = d_libMesh->elements_end();
            int N = 0;
            while ( pos != end ) {
                N++;
                ++pos;
            }
            n_ghost[i] = N - n_local[i];
        } else {
            // We are counting an intermediate element (not finished)
            n_local[i] = static_cast<size_t>(-1);
            n_global[i] = static_cast<size_t>(-1);
            n_ghost[i] = static_cast<size_t>(-1);
        }
    }
    // Compute the bounding box of the mesh
    d_box = std::vector<double>(PhysicalDim*2);
    for (int i=0; i<PhysicalDim; i++) {
        d_box[2*i+0] = 1e200;
        d_box[2*i+1] = -1e200;
    }
    ::Mesh::node_iterator node_pos = d_libMesh->nodes_begin();
    ::Mesh::node_iterator node_end = d_libMesh->nodes_end();
    while ( node_pos != node_end ) {
        ::Node* node = *node_pos;
        for (int i=0; i<PhysicalDim; i++) {
            double x = (*node)(i);
            if ( x < d_box[2*i+0] ) { d_box[2*i+0] = x; }
            if ( x > d_box[2*i+1] ) { d_box[2*i+1] = x; }
        }
        ++node_pos;
    }
    // Construct the element neighbor information
    d_libMesh->find_neighbors();
    // Construct the node neighbor information
    neighborNodeIDs = std::vector<unsigned int>(n_local[0],(unsigned int)-1);
    neighborNodes = std::vector< std::vector< ::Node* > >(n_local[0]);
    node_pos = d_libMesh->local_nodes_begin();
    node_end = d_libMesh->local_nodes_end();
    size_t i=0;
    while ( node_pos != node_end ) {
        ::Node *node = node_pos.operator*();
        neighborNodeIDs[i] = node->id();
        ++node_pos;
        i++;
    }
    AMP::Utilities::quicksort(neighborNodeIDs);
    ::Mesh::element_iterator elem_pos = d_libMesh->local_elements_begin();
    ::Mesh::element_iterator elem_end = d_libMesh->local_elements_end();
    std::vector< std::set<unsigned int> > tmpNeighborNodes(n_local[0]);
    int rank = d_comm.getRank();
    while ( elem_pos != elem_end ) {
        ::Elem *elem = elem_pos.operator*();
        for (i=0; i<elem->n_nodes(); i++) {
            ::Node *node = elem->get_node(i);
            int p_id = node->processor_id();
            if ( p_id==rank ) {
                int j = AMP::Utilities::findfirst(neighborNodeIDs,node->id());
                AMP_ASSERT(neighborNodeIDs[j]==node->id());
                for (unsigned int k=0; k<elem->n_nodes(); k++) {
                    if ( k==i )
                        continue;
                    ::Node *node2 = elem->get_node(k);
                    tmpNeighborNodes[j].insert(node2->id());
                }
            }
        }
        ++elem_pos;
    }
    for (i=0; i<neighborNodeIDs.size(); i++) {
        neighborNodes[i] = std::vector< ::Node* >(tmpNeighborNodes[i].size());
        int j = 0;
        for (std::set<unsigned int>::iterator it=tmpNeighborNodes[i].begin(); it!=tmpNeighborNodes[i].end(); it++) {
            neighborNodes[i][j] = d_libMesh->node_ptr(*it);
            j++;
        }
    }
    // Construct the boundary lists
    const std::set<short int> libmesh_bids = d_libMesh->boundary_info->get_boundary_ids();
    std::vector<short int> bids(libmesh_bids.begin(),libmesh_bids.end());
    Utilities::quicksort(bids);
    std::vector<short int> side_ids;
    std::vector<short int> node_ids;
    d_libMesh->boundary_info->build_side_boundary_ids(side_ids);
    d_libMesh->boundary_info->build_node_list_from_side_list();
    d_libMesh->boundary_info->build_node_boundary_ids(node_ids);
    for (int type2=0; type2<=(int)GeomDim; type2++) {
        GeomType type = (GeomType) type2;
        if ( type!=Vertex && type!=GeomDim ) {
            // Not all types are supported yet for iterators
            continue;
        }
        MeshIterator iterator = getIterator(type,0);
        for (size_t i=0; i<bids.size(); i++) {
            int id = (int) bids[i];
            // Count the number of elements on the given boundary
            MeshIterator curElem = iterator.begin();
            MeshIterator endElem = iterator.end();
            int N = 0;
            while ( curElem != endElem ) {
                if ( curElem->isOnBoundary(id) )
                    N++;
                ++curElem;
            }
            // Create the boundary list
            boost::shared_ptr<std::vector<MeshElement> > list( new std::vector<MeshElement>(N) );
            curElem = iterator.begin();
            endElem = iterator.end();
            N = 0;
            while ( curElem != endElem ) {
                if ( curElem->isOnBoundary(id) ) {
                    list->operator[](N) = *curElem;
                    N++;
                }
                ++curElem;
            }
            // Store the list
            std::pair<int,GeomType> mapid = std::pair<int,GeomType>(id,type);
            std::pair< std::pair<int,GeomType>, boost::shared_ptr<std::vector<MeshElement> > > entry(mapid,list);
            d_boundarySets.insert(entry);
        }
    }
}


/********************************************************
* Function to estimate the mesh size                    *
********************************************************/
size_t libMesh::estimateMeshSize( const MeshParameters::shared_ptr &params )
{
    boost::shared_ptr<AMP::Database> database = params->getDatabase();
    AMP_ASSERT(database.get()!=NULL);
    AMP_INSIST(database->keyExists("NumberOfElements"),"Key NumberOfElements must exist in database to estimate the mesh size");
    return (size_t) database->getInteger("NumberOfElements");
}


/********************************************************
* Return the number of elements                         *
********************************************************/
size_t libMesh::numLocalElements( const GeomType type ) const
{
    if ( n_local[type] == static_cast<size_t>(-1) )
        AMP_ERROR("numLocalElements is not implimented for this type");
    return n_local[type];
}
size_t libMesh::numGlobalElements( const GeomType type ) const
{
    if ( n_global[type] == static_cast<size_t>(-1) )
        AMP_ERROR("numLocalElements is not implimented for this type");
    return n_global[type];
}
size_t libMesh::numGhostElements( const GeomType type, int gcw ) const
{
    if ( gcw == 0 )
        return 0;
    if ( gcw > 1 )
        AMP_ERROR("Libmesh only supports a ghost cell width of 1");
    if ( n_ghost[type] == static_cast<size_t>(-1) )
        AMP_ERROR("numLocalElements is not implimented for this type");
    return n_ghost[type];
}


/********************************************************
* Return an iterator over the given geometric type      *
********************************************************/
MeshIterator libMesh::getIterator( const GeomType type, const int gcw )
{
    libMeshIterator iterator;
    if ( type==PhysicalDim ) {
        // This is a libMesh element
        if ( gcw==0 ) {
            ::Mesh::element_iterator begin = d_libMesh->local_elements_begin();
            ::Mesh::element_iterator end   = d_libMesh->local_elements_end();
            iterator = libMeshIterator( 1, this, gcw, (void*) &(begin), (void*) &(end), (void*) &(begin) );
        } else if ( gcw==1 ) {
            ::Mesh::element_iterator begin = d_libMesh->elements_begin();
            ::Mesh::element_iterator end   = d_libMesh->elements_end();
            iterator = libMeshIterator( 1, this, gcw, (void*) &(begin), (void*) &(end), (void*) &(begin) );
        } else {
            AMP_ERROR("Unsupported ghost cell width");
        }
    } else if ( type==Vertex ) {
        // This is a libMesh node
        if ( gcw==0 ) {
            ::Mesh::node_iterator begin = d_libMesh->local_nodes_begin();
            ::Mesh::node_iterator end   = d_libMesh->local_nodes_end();
            iterator = libMeshIterator( 0, this, gcw, (void*) &(begin), (void*) &(end), (void*) &(begin) );
        } else if ( gcw==1 ) {
            ::Mesh::node_iterator begin = d_libMesh->nodes_begin();
            ::Mesh::node_iterator end   = d_libMesh->nodes_end();
            iterator = libMeshIterator( 0, this, gcw, (void*) &(begin), (void*) &(end), (void*) &(begin) );
        } else {
            AMP_ERROR("Unsupported ghost cell width");
        }
    } else {
        AMP_ERROR("Not implimented for this type");
    }
    MeshIterator test = iterator;
    return test;
}


/********************************************************
* Return an iterator over the given boundary ids        *
* Note: we have not programmed this for ghosts yet      *
********************************************************/
std::vector<int> libMesh::getIDSets ( )
{
    const std::set<short int> libmesh_bids = d_libMesh->boundary_info->get_boundary_ids();
    std::vector<int> bids(libmesh_bids.size(),0);
    std::set<short int>::iterator it = libmesh_bids.begin();
    for (size_t i=0; i<bids.size(); i++) {
        bids[i] = *it;
        ++it;
    }
    return bids;
}
MeshIterator libMesh::getIDsetIterator ( const GeomType type, const int id, const int gcw )
{
    AMP_INSIST(gcw==0,"Iterator over ghost boundary elements is not supported yet");
    std::pair<int,GeomType> mapid = std::pair<int,GeomType>(id,type);
    std::map< std::pair<int,GeomType>, boost::shared_ptr<std::vector<MeshElement> > >::iterator it;
    boost::shared_ptr<std::vector<MeshElement> > list( new std::vector<MeshElement>() );
    it = d_boundarySets.find(mapid);
    if ( it != d_boundarySets.end() )
        list = it->second;
    return MultiVectorIterator( list, 0 );
}


/********************************************************
* Return pointers to the neighbor nodes give a node id  *
********************************************************/
std::vector< ::Node* > libMesh::getNeighborNodes( MeshElementID id )
{
    AMP_INSIST(id.type()==Vertex,"This function is for nodes");
    AMP_INSIST(id.meshID()==d_meshID,"Unknown mesh");
    AMP_INSIST(id.is_local(),"Only owned nodes can return their neighbor lists");
    int i = AMP::Utilities::findfirst(neighborNodeIDs,id.local_id());
    AMP_ASSERT(neighborNodeIDs[i]==id.local_id());
    return neighborNodes[i];
}


/********************************************************
* Displace a mesh by a scalar ammount                   *
********************************************************/
void libMesh::displaceMesh( std::vector<double> x_in )
{
    // Check x
    AMP_INSIST((short int)x_in.size()==PhysicalDim,"Displacement vector size should match PhysicalDim");
    std::vector<double> x = x_in;
    d_comm.minReduce(&x[0],x.size());
    for (size_t i=0; i<x.size(); i++)
        AMP_INSIST(fabs(x[i]-x_in[i])<1e-12,"x does not match on all processors");
    // Move the mesh
    ::Mesh::node_iterator cur = d_libMesh->nodes_begin();
    ::Mesh::node_iterator end = d_libMesh->nodes_end();
    while ( cur != end ) {
        ::Node  *d_Node = *cur;
        for (size_t i=0; i<x.size(); i++)
            (*d_Node)(i) += x[i];
        ++cur;
    }
    // Update the bounding box
    for (int i=0; i<PhysicalDim; i++) {
        d_box[2*i+0] += x[i];
        d_box[2*i+1] += x[i];
    }
}


} // Mesh namespace
} // AMP namespace

