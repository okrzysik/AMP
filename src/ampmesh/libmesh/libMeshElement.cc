#include "ampmesh/libmesh/libMeshElement.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {


/********************************************************
* Constructors                                          *
********************************************************/
libMeshElement::libMeshElement()
{
}
libMeshElement::libMeshElement(int dim, int type, void* element)
{
    d_dim = dim;
    d_type = type;
    ptr_element = element;
    if ( type==0 ) {
        d_elementType = Vertex;
        ::Node* node = (::Node*) ptr_element;
        d_globalID = node->id();
    } else if ( type==1 ) {
        d_elementType = (GeomType) dim;
        ::Elem* elem = (::Elem*) ptr_element;
        d_globalID = elem->id();
    } else {
        AMP_ERROR("libMesh does not this element type");
    }
}


/********************************************************
* De-constructor                                        *
********************************************************/
libMeshElement::~libMeshElement()
{
}





} // Mesh namespace
} // AMP namespace

