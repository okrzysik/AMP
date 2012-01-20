#include "ampmesh/libmesh/libMeshElement.h"
#include "utils/Utilities.h"

#include "boundary_info.h"

namespace AMP {
namespace Mesh {


// Create a unique id for this class
static unsigned int libMeshElementTypeID = TYPE_HASH(libMeshElement);


/********************************************************
* Constructors                                          *
********************************************************/
libMeshElement::libMeshElement()
{
    typeID = libMeshElementTypeID;
    element = NULL;
    d_dim = -1;
    d_elementType = null;
    d_globalID = MeshElementID();
}
libMeshElement::libMeshElement(int dim, GeomType type, void* libmesh_element, 
    unsigned int rank, MeshID meshID, libMesh* mesh)
{
    AMP_ASSERT(libmesh_element!=NULL);
    typeID = libMeshElementTypeID;
    element = NULL;
    d_elementType = type;
    d_dim = dim;
    d_rank = rank;
    d_mesh = mesh;
    d_meshID = meshID;
    ptr_element = libmesh_element;
    unsigned int local_id = (unsigned int)-1;
    unsigned int owner_rank = (unsigned int)-1;
    bool is_local=false;
    if ( d_elementType==Vertex ) {
        d_elementType = Vertex;
        ::Node* node = (::Node*) ptr_element;
        local_id = node->id();
        owner_rank = node->processor_id();
        is_local = owner_rank==d_rank;
    } else {
        d_elementType = (GeomType) dim;
        ::Elem* elem = (::Elem*) ptr_element;
        AMP_ASSERT(elem->n_neighbors()<100);
        local_id = elem->id();
        owner_rank = elem->processor_id();
        is_local = owner_rank==d_rank;
    }
    d_globalID = MeshElementID(is_local,d_elementType,local_id,owner_rank,meshID);
}
libMeshElement::libMeshElement(int dim, GeomType type, boost::shared_ptr< ::Elem > libmesh_element, 
    unsigned int rank, MeshID meshID, libMesh* mesh)
{
    AMP_ASSERT(libmesh_element.get()!=NULL);
    typeID = libMeshElementTypeID;
    element = NULL;
    d_elementType = type;
    d_dim = dim;
    d_rank = rank;
    d_mesh = mesh;
    d_meshID = meshID;
    ptr2 = libmesh_element;
    ptr_element = libmesh_element.get();
    unsigned int local_id = (unsigned int)-1;
    unsigned int owner_rank = (unsigned int)-1;
    bool is_local=false;
    if ( d_elementType==Vertex ) {
        ::Node* node = (::Node*) ptr_element;
        local_id = node->id();
        owner_rank = node->processor_id();
        is_local = owner_rank==d_rank;
    } else {
        ::Elem* elem = (::Elem*) ptr_element;
        local_id = elem->id();
        owner_rank = elem->processor_id();
        is_local = owner_rank==d_rank;
    }
    d_globalID = MeshElementID(is_local,d_elementType,local_id,owner_rank,meshID);
}
libMeshElement::libMeshElement(const libMeshElement& rhs)
{
    typeID = libMeshElementTypeID;
    element = NULL;
    d_elementType = rhs.d_elementType;
    d_globalID = rhs.d_globalID;
    d_dim = rhs.d_dim;
    ptr_element = rhs.ptr_element;
    ptr2 = rhs.ptr2;
    d_rank = rhs.d_rank;
    d_mesh = rhs.d_mesh;
    d_meshID = rhs.d_meshID;
}
libMeshElement& libMeshElement::operator=(const libMeshElement& rhs)
{
    if (this == &rhs) // protect against invalid self-assignment
        return *this;
    this->typeID = libMeshElementTypeID;
    this->element = NULL;
    this->d_elementType = rhs.d_elementType;
    this->d_globalID = rhs.d_globalID;
    this->d_dim = rhs.d_dim;
    this->ptr_element = rhs.ptr_element;
    this->ptr2 = rhs.ptr2;
    this->d_rank = rhs.d_rank;
    this->d_mesh = rhs.d_mesh;
    this->d_meshID = rhs.d_meshID;
    return *this;
}


/****************************************************************
* De-constructor                                                *
****************************************************************/
libMeshElement::~libMeshElement()
{
    element = NULL;
}


/****************************************************************
* Function to clone the element                                 *
****************************************************************/
MeshElement* libMeshElement::clone() const
{
    return new libMeshElement(*this);
}


/****************************************************************
* Function to get the elements composing the current element    *
****************************************************************/
std::vector<MeshElement> libMeshElement::getElements(const GeomType type) const
{
    AMP_INSIST(type<=d_elementType,"sub-elements must be of a smaller or equivalent type");
    std::vector<MeshElement> children(0);
    ::Elem* elem = (::Elem*) ptr_element;
    if ( d_elementType==Vertex ) {
        // A vertex does not have children, return itself
        if ( type!=Vertex )
            AMP_ERROR("A vertex is the base element and cannot have and sub-elements");
        children.resize(1);
        children[0] = *this;
    } else if ( type==d_elementType ) {
        // Return the children of the current element
        if ( elem->has_children() ) {
            children.resize(elem->n_children());
            for (unsigned int i=0; i<children.size(); i++)
                children[i] = libMeshElement( d_dim, type, (void*)elem->child(i), d_rank, d_meshID, d_mesh );
        } else {
            children.resize(1);
            children[0] = *this;
        }
    } else if ( type==Vertex ) {
        // Return the nodes of the current element
        children.resize(elem->n_nodes());
        for (unsigned int i=0; i<children.size(); i++)
            children[i] = libMeshElement( d_dim, type, (void*)elem->get_node(i), d_rank, d_meshID, d_mesh );
    } else if ( ((int) d_elementType - (int) type) == 1 ) {
        // Return the faces of the current element
        children.resize(elem->n_faces());
        std::vector< MeshElement::shared_ptr > neighbors = this->getNeighbors();
        for (unsigned int i=0; i<children.size(); i++) {
            // We need to build a valid element
            ::AutoPtr< ::Elem > tmp = elem->build_side(i,false);
            boost::shared_ptr< ::Elem > element( tmp.release() );
            // Create the id and owning processor for the element
            MeshElementID global_id = d_globalID;
            if ( neighbors[i].get() != NULL ) {
                if ( neighbors[i]->globalID() < global_id )
                    global_id = neighbors[i]->globalID();
            }
            AMP_ASSERT(children.size()<32);
            unsigned int id = (global_id.local_id()<<5) + i;
            element->set_id() = id;
            element->processor_id() = global_id.owner_rank();
            // Create the libMeshElement
            children[i] = libMeshElement( d_dim, type, element, d_rank, d_meshID, d_mesh );
        }
    } else if  ( ((int) d_elementType - (int) type) == 2 ) {
        // Return the edges of the current element
        children.resize(elem->n_edges());
        for (unsigned int i=0; i<children.size(); i++) {
            // We need to build a valid element
            ::AutoPtr< ::Elem > tmp = elem->build_edge(i);
            boost::shared_ptr< ::Elem > element( tmp.release() );
            // Create the id and owning processor for the element
            if ( type == Edge ) {
                // Use the two nodes to construct the id and owning processor
                AMP_ASSERT(element->n_nodes()==2);
                ::Node *node1 = element->get_node(0);
                ::Node *node2 = element->get_node(1);
                unsigned int node1_id =  node1->id();
                unsigned int node2_id =  node2->id();
                if ( node1_id < node2_id ) {
                    element->processor_id() = node1->processor_id();
                    element->set_id() = (node1_id<<15) + node2_id;
                } else {
                    element->processor_id() = node2->processor_id();
                    element->set_id() = (node2_id<<15) + node1_id;
                }
            } else {
                // Not what?
                AMP_ERROR("Internal error in libMeshElement::getElements");
            }
            // Create the libMeshElement
            children[i] = libMeshElement( d_dim, type, element, d_rank, d_meshID, d_mesh );
        }

    }
    return children;
}


/****************************************************************
* Function to get the neighboring elements                      *
****************************************************************/
std::vector< MeshElement::shared_ptr > libMeshElement::getNeighbors() const
{
    std::vector< MeshElement::shared_ptr > neighbors(0);
    if ( d_elementType==Vertex ) {
        // Return the neighbors of the current node
        std::vector< ::Node* > neighbor_nodes = d_mesh->getNeighborNodes( d_globalID );
        neighbors.resize(neighbor_nodes.size(),MeshElement::shared_ptr());
        for (size_t i=0; i<neighbor_nodes.size(); i++) {
            // There are no NULL neighbors
            boost::shared_ptr<libMeshElement> neighbor(new libMeshElement( d_dim, Vertex, (void*)neighbor_nodes[i], d_rank, d_meshID, d_mesh ) );
            neighbors[i] = neighbor;
        }
    } else {
        // Return the neighbors of the current element
        ::Elem* elem = (::Elem*) ptr_element;
        //if ( elem->n_neighbors()==0 )
        //    AMP_ERROR("Element has not neighbors, this could indicate a problem with the mesh");
    	neighbors.resize(elem->n_neighbors());
        for (size_t i=0; i<neighbors.size(); i++) {
            void *neighbor_elem = (void*) elem->neighbor(i);
            boost::shared_ptr<libMeshElement> neighbor;
            if ( neighbor_elem != NULL )
                neighbor = boost::shared_ptr<libMeshElement>(new libMeshElement( d_dim, d_elementType, neighbor_elem, d_rank, d_meshID, d_mesh ) );
            neighbors[i] = neighbor;
        }
    }
    return neighbors;
}


/****************************************************************
* Functions to get basic element properties                     *
****************************************************************/
double libMeshElement::volume() const
{
    if ( d_elementType == Vertex )
        AMP_ERROR("volume is is not defined Nodes");
    ::Elem* elem = (::Elem*) ptr_element;
    return elem->volume();
}
std::vector<double> libMeshElement::coord() const
{
    if ( d_elementType != Vertex )
        AMP_ERROR("coord is only defined for Nodes");
    ::Node* node = (::Node*) ptr_element;
    std::vector<double> x(d_dim,0.0);
    for (int i=0; i<d_dim; i++)
        x[i] = (*node)(i);
    return x;
}
bool libMeshElement::isOnSurface() const
{
    MeshElement search = MeshElement(*this);
    std::vector<MeshElement> &data = *(d_mesh->d_surfaceSets[d_elementType]);
    size_t index = AMP::Utilities::findfirst( data, search );
    if ( index < data.size() ) {
        if ( d_mesh->d_surfaceSets[d_elementType]->operator[](index).globalID() == d_globalID )
            return true;
    }
    return false;
}
bool libMeshElement::isOnBoundary(int id) const
{
    bool on_boundary = false;
    boost::shared_ptr< ::Mesh> d_libMesh = d_mesh->getlibMesh();
    if ( d_elementType==Vertex ) {
        // Entity is a libmesh node
        ::Node* node = (::Node*) ptr_element;
        std::vector< short int > bids = d_libMesh->boundary_info->boundary_ids(node);
        for (size_t i=0; i<bids.size(); i++) {
            if ( bids[i]==id )
                on_boundary = true;
        }
    } else if ( (int)d_elementType==d_dim ) {
        // Entity is a libmesh node
        ::Elem* elem = (::Elem*) ptr_element;
        unsigned int side = d_libMesh->boundary_info->side_with_boundary_id(elem,id);
        if ( side != static_cast<unsigned int>(-1) )
            on_boundary = true;
    } else  {
        // All other entities are on the boundary iff all of their verticies are on the surface
        std::vector<MeshElement> nodes = this->getElements(Vertex);
        on_boundary = true;
        for (size_t i=0; i<nodes.size(); i++) 
            on_boundary = on_boundary && nodes[i].isOnBoundary(id);
    }
    return on_boundary;
}

} // Mesh namespace
} // AMP namespace

