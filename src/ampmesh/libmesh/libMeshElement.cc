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
    d_globalID = 0;
}
libMeshElement::libMeshElement(int dim, GeomType type, void* libmesh_element)
{
    typeID = libMeshElementTypeID;
    element = NULL;
    d_elementType = type;
    d_dim = dim;
    ptr_element = libmesh_element;
    if ( d_elementType==Vertex ) {
        d_elementType = Vertex;
        ::Node* node = (::Node*) ptr_element;
        d_globalID = node->id();
    } else {
        d_elementType = (GeomType) dim;
        ::Elem* elem = (::Elem*) ptr_element;
        d_globalID = elem->id();
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
std::vector<MeshElement> libMeshElement::getElements(GeomType &type)
{
    if ( d_elementType==Vertex )
        AMP_ERROR("A vertex is the base element and cannot have and sub-elements");
    if ( type >= d_elementType )
        AMP_ERROR("type must be <= the GeomType of the current element");
    std::vector<MeshElement> children(0);
    ::Elem* elem = (::Elem*) ptr_element;
    if ( type==d_elementType ) {
        // Return the children of the current element
        children.resize(elem->n_children());
        for (unsigned int i=0; i<children.size(); i++)
            children[i] = libMeshElement( d_dim, type, (void*)elem->child(i) );
    } else if ( type==Vertex ) {
        // Return the nodes of the current element
        children.resize(elem->n_nodes());
        for (unsigned int i=0; i<children.size(); i++)
            children[i] = libMeshElement( d_dim, type, (void*)elem->get_node(i) );
    } else if ( type==Edge ) {
        // Return the edges of the current element
        AMP_ERROR("unfinished");
        //children.resize(elem->n_edges());
        //for (unsigned int i=0; i<children.size(); i++)
        //    children[i] = libMeshElement( d_dim, type, (void*)elem->build_edge(i) );
    } else if ( type==Face ) {
        // Return the edges of the current element
        AMP_ERROR("unfinished");
        //children.resize(elem->n_faces());
        //for (unsigned int i=0; i<children.size(); i++)
        //    //children[i] = libMeshElement( d_dim, type, (void*)elem->build_face(i) );
    }
    return children;
}
std::vector<MeshElement> libMeshElement::getNeighbors()
{
    AMP_ERROR("getNeighbors is not finished");
    return std::vector<MeshElement>(0);
}
double libMeshElement::volume()
{
    if ( d_elementType == Vertex )
        AMP_ERROR("volume is is not defined Nodes");
    ::Elem* elem = (::Elem*) ptr_element;
    return elem->volume();
}
std::vector<double> libMeshElement::coord()
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

