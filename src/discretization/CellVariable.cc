#include "CellVariable.h"
#include "ampmesh/MeshElement.h"

namespace AMP {
namespace Discretization {


/********************************************************
* Constructors                                          *
********************************************************/
CellVariable::CellVariable( int DOFsPerCell, const std::string &name ):
    AMP::LinearAlgebra::Variable::Variable(name)
{
    d_VariableName = name;
    d_DOFsPerCell = DOFsPerCell;
}


/********************************************************
* Return number of DOFs per node                        *
********************************************************/
size_t CellVariable::DOFsPerObject() const
{
    return d_DOFsPerCell;
}


/********************************************************
* Return the variable ID type                           *
********************************************************/
size_t CellVariable::variableID() const
{
    return (size_t) AMP::Mesh::Vertex;
}


}
}

