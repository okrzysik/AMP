#include "ampmesh/libmesh/libMeshElement.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {

// Create a unique id for this class
static int libMeshElementTypeID = 2;


/********************************************************
* Constructors                                          *
********************************************************/
libMeshElement::libMeshElement()
{
    typeID = libMeshElementTypeID;
    element = NULL;
}
libMeshElement::libMeshElement(int dim, int type, void* libmesh_element)
{
    typeID = libMeshElementTypeID;
    d_dim = dim;
    d_type = type;
    element = NULL;
    ptr_element = libmesh_element;
    if ( d_type==0 ) {
        d_elementType = Vertex;
        ::Node* node = (::Node*) ptr_element;
        d_globalID = node->id();
    } else if ( d_type==1 ) {
        d_elementType = (GeomType) dim;
        ::Elem* elem = (::Elem*) ptr_element;
        d_globalID = elem->id();
    } else {
        AMP_ERROR("libMesh does not define this element type");
    }
}
libMeshElement::libMeshElement(const libMeshElement& rhs)
{
    typeID = libMeshElementTypeID;
    element = NULL;
    d_dim = rhs.d_dim;
    d_type = rhs.d_type;
    ptr_element = rhs.ptr_element;
    d_elementType = rhs.d_elementType;
    d_globalID = rhs.d_globalID;
}
libMeshElement& libMeshElement::operator=(const libMeshElement& rhs)
{
    if (this == &rhs) // protect against invalid self-assignment
        return *this;
    this->typeID = libMeshElementTypeID;
    this->element = NULL;
    this->d_dim = rhs.d_dim;
    this->d_type = rhs.d_type;
    this->ptr_element = rhs.ptr_element;
    this->d_elementType = rhs.d_elementType;
    this->d_globalID = rhs.d_globalID;
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
    return std::vector<MeshElement>(0);
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

