#include "NodalVariable.h"
#include "ampmesh/MeshElement.h"

namespace AMP {
namespace Discretization {


/********************************************************
* Constructors                                          *
********************************************************/
NodalVariable::NodalVariable( int DOFsPerNode, const std::string &name ):
    AMP::LinearAlgebra::Variable::Variable(name)
{
    d_VariableName = name;
    d_DOFsPerNode = DOFsPerNode;
}


/********************************************************
* Return number of DOFs per node                        *
********************************************************/
size_t NodalVariable::DOFsPerObject() const
{
    return d_DOFsPerNode;
}


/********************************************************
* Return the variable ID type                           *
********************************************************/
size_t NodalVariable::variableID() const
{
    return (size_t) AMP::Mesh::Vertex;
}


}
}

