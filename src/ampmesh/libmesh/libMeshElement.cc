#include "ampmesh/libmesh/libMeshElement.h"
#include "utils/Utilities.h"

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
    unsigned int rank, size_t meshID, libMesh* mesh)
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
    d_globalID = MeshElementID();
    d_globalID.type = type;
    d_globalID.meshID = meshID;
    if ( d_elementType==Vertex ) {
        d_elementType = Vertex;
        ::Node* node = (::Node*) ptr_element;
        d_globalID.local_id = node->id();
        d_globalID.owner_rank = node->processor_id();
        d_globalID.is_local = d_globalID.owner_rank==d_rank;
    } else {
        d_elementType = (GeomType) dim;
        ::Elem* elem = (::Elem*) ptr_element;
        d_globalID.local_id = elem->id();
        d_globalID.owner_rank = elem->processor_id();
        d_globalID.is_local = d_globalID.owner_rank==d_rank;
    }
}
libMeshElement::libMeshElement(const libMeshElement& rhs)
{
    typeID = libMeshElementTypeID;
    element = NULL;
    d_elementType = rhs.d_elementType;
    d_globalID = rhs.d_globalID;
    d_dim = rhs.d_dim;
    ptr_element = rhs.ptr_element;
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
    this->d_rank = rhs.d_rank;
    this->d_mesh = rhs.d_mesh;
    this->d_meshID = rhs.d_meshID;
    return *this;
}


/********************************************************
* De-constructor                                        *
********************************************************/
libMeshElement::~libMeshElement()
{
    element = NULL;
}


/********************************************************
* Function to clone the element                         *
********************************************************/
MeshElement* libMeshElement::clone() const
{
    return new libMeshElement(*this);
}


/********************************************************
* Functions that aren't implimented yet                 *
********************************************************/
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
    } else if ( type==Edge ) {
        // Return the edges of the current element
        AMP_ERROR("unfinished");
        //children.resize(elem->n_edges());
        //for (unsigned int i=0; i<children.size(); i++)
        //    children[i] = libMeshElement( d_dim, type, (void*)elem->build_edge(i), d_rank, d_meshID, d_mesh );
    } else if ( type==Face ) {
        // Return the edges of the current element
        AMP_ERROR("unfinished");
        //children.resize(elem->n_faces());
        //for (unsigned int i=0; i<children.size(); i++)
        //    //children[i] = libMeshElement( d_dim, type, (void*)elem->build_face(i), d_rank, d_meshID, d_mesh );
    }
    return children;
}
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
    } else if ( (int)d_elementType==d_dim ) {
        // Return the neighbors of the current element
        ::Elem* elem = (::Elem*) ptr_element;
    	neighbors.resize(elem->n_neighbors());
        for (size_t i=0; i<neighbors.size(); i++) {
            void *neighbor_elem = (void*) elem->neighbor(i);
            boost::shared_ptr<libMeshElement> neighbor;
            if ( neighbor_elem != NULL )
                neighbor = boost::shared_ptr<libMeshElement>(new libMeshElement( d_dim, d_elementType, neighbor_elem, d_rank, d_meshID, d_mesh ) );
            neighbors[i] = neighbor;
        }
    } else  {
        AMP_ERROR("unfinished");
    }
    return neighbors;
}
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



} // Mesh namespace
} // AMP namespace

