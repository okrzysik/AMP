#include "ampmesh/MeshElement.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Mesh {


/********************************************************
* Constructors                                          *
********************************************************/
MeshElement::MeshElement()
{
}


/********************************************************
* De-constructor                                        *
********************************************************/
MeshElement::~MeshElement()
{
}


/********************************************************
* Functions that aren't implimented for the base class  *
********************************************************/
std::vector<MeshElement> MeshElement::getElements(GeomType &type)
{
    AMP_ERROR("Not Implimented Yet");
    return std::vector<MeshElement>(0);
}
std::vector<MeshElement> MeshElement::getNeighbors()
{
    AMP_ERROR("Not Implimented Yet");
    return std::vector<MeshElement>(0);
}
double MeshElement::getVolume()
{
    AMP_ERROR("Not Implimented Yet");
    return 0.0;
}
void MeshElement::getCoord(std::vector<double> &coords)
{
    AMP_ERROR("Not Implimented Yet");
}


} // Mesh namespace
} // AMP namespace

