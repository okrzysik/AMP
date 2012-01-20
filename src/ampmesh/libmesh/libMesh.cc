#include "ampmesh/libmesh/initializeLibMesh.h"
#include "ampmesh/libmesh/libMesh.h"
#include "ampmesh/libmesh/libMeshElement.h"
#include "ampmesh/libmesh/libMeshIterator.h"
#include "ampmesh/MeshElementVectorIterator.h"
#include "ampmesh/MultiIterator.h"
#include "utils/MemoryDatabase.h"
#include "utils/AMPManager.h"
#include "utils/Utilities.h"

#ifdef USE_AMP_VECTORS
    #include "vectors/Vector.h"
    #include "vectors/Variable.h"
    #include "vectors/VectorBuilder.h"
#endif
#ifdef USE_AMP_DISCRETIZATION
    #include "discretization/DOF_Manager.h"
    #include "discretization/simpleDOF_Manager.h"
#endif


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
    // First we need to destroy the elements, surface sets, and boundary sets
    d_localElements.clear();
    d_ghostElements.clear();
    d_surfaceSets.clear();
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
    // Construct the list of elements of type side or edge
    for (int i=0; i<=(int)GeomDim; i++) {
        GeomType type = (GeomType) i;
        if ( type==Vertex || type==GeomDim )
            continue;
        // Get a unique list of all elements of the desired type
        // Note: we will include ghost elements of the desired type, but only
        // those elements that are part of an owned element of type GeomDim.
        // This is required because of the way we create the id of the new element
        std::set<MeshElement> element_list;
        MeshIterator it = getIterator( GeomDim, 0 );
        for (size_t j=0; j<it.size(); j++) {
            std::vector<MeshElement> tmp = it->getElements(type);
            for (size_t k=0; k<tmp.size(); k++)
                element_list.insert(tmp[k]);
            ++it;
        }
        // Split the new elements into the local and ghost lists
        size_t N_local = 0;
        size_t N_ghost = 0;
        for (std::set<MeshElement>::iterator it2 = element_list.begin(); it2!=element_list.end(); it2++) {
            MeshElementID id = it2->globalID();
            if ( id.is_local() )
                N_local++;
            else
                N_ghost++;
        }
        size_t N_global = d_comm.sumReduce(N_local);
        AMP_ASSERT(N_global>=n_global[GeomDim]);
        boost::shared_ptr<std::vector<MeshElement> >  local_elements( new std::vector<MeshElement>(N_local) );
        boost::shared_ptr<std::vector<MeshElement> >  ghost_elements( new std::vector<MeshElement>(N_ghost) );
        N_local = 0;
        N_ghost = 0;
        for (std::set<MeshElement>::iterator it2 = element_list.begin(); it2!=element_list.end(); it2++) {
            MeshElementID id = it2->globalID();
            if ( id.is_local() ) {
                local_elements->operator[](N_local) = *it2;
                N_local++;
            } else {
                ghost_elements->operator[](N_ghost) = *it2;
                N_ghost++;
            }
        }
        std::pair< GeomType, boost::shared_ptr<std::vector<MeshElement> > > local_pair( type, local_elements );
        std::pair< GeomType, boost::shared_ptr<std::vector<MeshElement> > > ghost_pair( type, ghost_elements );
        d_localElements.insert( local_pair );
        d_ghostElements.insert( ghost_pair );
    }
    // Construct the boundary elements for Node and Elem
    d_surfaceSets = std::vector< boost::shared_ptr<std::vector<MeshElement> > >((int)GeomDim+1);
    elem_pos = d_libMesh->local_elements_begin();
    elem_end = d_libMesh->local_elements_end();
    std::set< ::Elem* > boundaryElements;
    std::set< ::Node* > boundaryNodes;
    while ( elem_pos != elem_end ) {
        ::Elem* element = *elem_pos;
        if ( element->on_boundary() ) {
            boundaryElements.insert(element);
            for (unsigned int i=0; i<element->n_sides(); i++) {
                if ( element->neighbor(i) == NULL ) {
                    ::AutoPtr< ::Elem > side = element->build_side(i);
                    for (unsigned int j=0; j<side->n_nodes(); j++)
                        boundaryNodes.insert(side->get_node(j));
                }
            }
        }
        elem_pos++;
    }
    d_surfaceSets[GeomDim] = boost::shared_ptr<std::vector<MeshElement> >( new std::vector<MeshElement>(boundaryElements.size()) );
    std::set< ::Elem* >::iterator elem_iterator= boundaryElements.begin();
    for (size_t i=0; i<boundaryElements.size(); i++) {
        (*d_surfaceSets[GeomDim])[i] = libMeshElement(PhysicalDim, GeomDim, (void*) *elem_iterator, d_comm.getRank(), d_meshID, this );
        elem_iterator++;
    }
    AMP::Utilities::quicksort(*d_surfaceSets[GeomDim]);
    d_surfaceSets[Vertex] = boost::shared_ptr<std::vector<MeshElement> >( new std::vector<MeshElement>(boundaryNodes.size()) );
    std::set< ::Node* >::iterator node_iterator = boundaryNodes.begin();
    for (size_t i=0; i<boundaryNodes.size(); i++) {
        (*d_surfaceSets[Vertex])[i] = libMeshElement(PhysicalDim, Vertex, (void*) *node_iterator, d_comm.getRank(), d_meshID, this );
        node_iterator++;
    }
    AMP::Utilities::quicksort(*d_surfaceSets[Vertex]);
    // Construct the boundary elements for all other types
    // An face or edge is on the boundary if all of its nodes are on the surface
    size_t element_surface_global_size = d_comm.sumReduce(d_surfaceSets[GeomDim]->size());
    for (int type2=1; type2<(int)GeomDim; type2++) {
        GeomType type = (GeomType) type2;
        std::set<MeshElement> list;
        MeshIterator it = getIterator( type, 0 );
        for (size_t i=0; i<it.size(); i++) {
            std::vector<MeshElement> nodes = it->getElements(Vertex);
            AMP_ASSERT(nodes.size()>0);
            bool on_boundary = true;
            for (size_t j=0; j<nodes.size(); j++) {
                if ( !nodes[j].isOnSurface() )
                    on_boundary = false;
            }
            if ( on_boundary ) 
                list.insert(*it);
            ++it;
        }
        d_surfaceSets[type2] = boost::shared_ptr<std::vector<MeshElement> >( 
            new std::vector<MeshElement>(list.begin(),list.end()) );
        AMP::Utilities::quicksort(*d_surfaceSets[type2]);
        size_t local_size = d_surfaceSets[type2]->size();
        size_t global_size = d_comm.sumReduce(local_size);
        AMP_ASSERT(global_size>=element_surface_global_size);
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
    // Finish counting the elements 
    for (int i=0; i<=(int)GeomDim; i++) {
        GeomType type = (GeomType) i;
        if ( type==Vertex || type==GeomDim )
            continue;
        MeshIterator it = getIterator( type, 0 );
        n_local[i] = it.size();
        n_global[i] = d_comm.sumReduce(n_local[i]);
        it = getIterator( type, 1 );
        n_ghost[i] = it.size() - n_local[i];
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
        // All other types require a pre-constructed list
        std::map< GeomType, boost::shared_ptr<std::vector<MeshElement> > >::iterator it1, it2;
        if ( gcw==0 ) {
            it1 = d_localElements.find( type );
            if ( it1==d_localElements.end() )
                AMP_ERROR("Internal error in libMesh::getIterator");
            return MultiVectorIterator( it1->second, 0 );
        } else if ( gcw==1 ) {
            it1 = d_localElements.find( type );
            it2 = d_ghostElements.find( type );
            if ( it1==d_localElements.end() || it2==d_ghostElements.end() )
                AMP_ERROR("Internal error in libMesh::getIterator");
            std::vector<boost::shared_ptr<MeshIterator> > iterators(2);
            iterators[0] = boost::shared_ptr<MeshIterator>( new MultiVectorIterator( it1->second, 0 ) );
            iterators[1] = boost::shared_ptr<MeshIterator>( new MultiVectorIterator( it2->second, 0 ) );
            return MultiIterator( iterators, 0 );
        } else {
            AMP_ERROR("Unsupported ghost cell width");
        }
    }
    MeshIterator test = iterator;
    return test;
}


/********************************************************
* Return an iterator over the given boundary ids        *
* Note: we have not programmed this for ghosts yet      *
********************************************************/
MeshIterator libMesh::getSurfaceIterator ( const GeomType type, const int gcw )
{
    AMP_INSIST(gcw==0,"Iterator over ghost boundary elements is not supported yet");
    AMP_ASSERT( type>=0 && type<=GeomDim );
    boost::shared_ptr<std::vector<MeshElement> > list = d_surfaceSets[type];
    if ( list.get() == NULL )
        AMP_ERROR("Surface iterator over the given geometry type is not supported");
    return MultiVectorIterator( list, 0 );
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
* Return the position vector                            *
********************************************************/
#ifdef USE_AMP_VECTORS
AMP::LinearAlgebra::Vector::shared_ptr  libMesh::getPositionVector( std::string name, const int gcw )
{
    #ifdef USE_AMP_DISCRETIZATION
        AMP::Discretization::DOFManager::shared_ptr DOFs = 
            AMP::Discretization::simpleDOFManager::create( 
            shared_from_this(), AMP::Mesh::Vertex, gcw, PhysicalDim, false );
        AMP::LinearAlgebra::Variable::shared_ptr nodalVariable( new AMP::LinearAlgebra::Variable(name) );
        AMP::LinearAlgebra::Vector::shared_ptr position = AMP::LinearAlgebra::createVector( DOFs, nodalVariable, false );
        std::vector<size_t> dofs(PhysicalDim);
        AMP::Mesh::MeshIterator cur = DOFs->getIterator();
        AMP::Mesh::MeshIterator end = cur.end();
        while ( cur != end ) {
            AMP::Mesh::MeshElementID id = cur->globalID();
            std::vector<double> coord = cur->coord();
            DOFs->getDOFs( id, dofs );
            position->setValuesByGlobalID( dofs.size(), &dofs[0], &coord[0] );
            ++cur;
        }
        return position;
    #else
        AMP_ERROR("getPositionVector requires DISCRETIZATION");
        return AMP::LinearAlgebra::Vector::shared_ptr();
    #endif
}
#endif


/********************************************************
* Displace a mesh                                       *
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
#ifdef USE_AMP_VECTORS
void libMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr x )
{
    #ifdef USE_AMP_DISCRETIZATION
        // Create the position vector with the necessary ghost nodes
        AMP::Discretization::DOFManager::shared_ptr DOFs = 
        AMP::Discretization::simpleDOFManager::create( 
            shared_from_this(), getIterator(AMP::Mesh::Vertex,1), getIterator(AMP::Mesh::Vertex,0), PhysicalDim );
        AMP::LinearAlgebra::Variable::shared_ptr nodalVariable( new AMP::LinearAlgebra::Variable( "tmp_pos" ) );
        AMP::LinearAlgebra::Vector::shared_ptr displacement = AMP::LinearAlgebra::createVector( DOFs, nodalVariable, false );
        std::vector<size_t> dofs1(PhysicalDim);
        std::vector<size_t> dofs2(PhysicalDim);
        AMP::Mesh::MeshIterator cur = getIterator(AMP::Mesh::Vertex,0);
        AMP::Mesh::MeshIterator end = cur.end();
        AMP::Discretization::DOFManager::shared_ptr DOFx = x->getDOFManager();
        std::vector<double> data(PhysicalDim);
        while ( cur != end ) {
            AMP::Mesh::MeshElementID id = cur->globalID();
            DOFx->getDOFs( id, dofs1 );
            DOFs->getDOFs( id, dofs2 );
            x->getValuesByGlobalID( PhysicalDim, &dofs1[0], &data[0] );
            displacement->setValuesByGlobalID( PhysicalDim, &dofs2[0], &data[0] );
            ++cur;
        }
        displacement->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        // Move all nodes (including the ghost nodes
        ::Mesh::node_iterator node_cur = d_libMesh->nodes_begin();
        ::Mesh::node_iterator node_end = d_libMesh->nodes_end();
        int rank = d_comm.getRank();
        while ( node_cur != node_end ) {
            ::Node* node = *node_cur;
            // Create the element id
            unsigned int owner_rank = node->processor_id();
            unsigned int local_id = node->id();
            bool is_local = (int)owner_rank==rank;
            AMP::Mesh::MeshElementID id( is_local ,AMP::Mesh::Vertex, local_id, owner_rank, d_meshID);
            // Get the position of the point
            DOFs->getDOFs( id, dofs2 );
            displacement->getValuesByGlobalID( PhysicalDim, &dofs2[0], &data[0] );
            // Move the point
            for (int i=0; i<PhysicalDim; i++)
                (*node)(i) += data[i];
            ++node_cur;
        }
    #else
        AMP_ERROR("getPositionVector requires DISCRETIZATION");
    #endif
}
#endif


} // Mesh namespace
} // AMP namespace

