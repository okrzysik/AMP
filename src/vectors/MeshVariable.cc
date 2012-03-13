#ifdef USE_AMP_MESH
#include "MeshVariable.h"

namespace AMP {
namespace LinearAlgebra {


MeshVariable::MeshVariable ( const std::string &name, AMP::Mesh::Mesh::shared_ptr mesh, bool useMeshComm ):
    SubsetVariable ( name )
{
    AMP_ASSERT(mesh.get()!=NULL);
    d_mesh = mesh;
    d_useMeshComm = useMeshComm;
}


AMP::Discretization::DOFManager::shared_ptr  MeshVariable::getSubsetDOF( AMP::Discretization::DOFManager::shared_ptr parentDOF )
{
    return parentDOF->subset( d_mesh, d_useMeshComm );
}


MeshIteratorVariable::MeshIteratorVariable ( const std::string &name, const AMP::Mesh::MeshIterator &iterator, const AMP_MPI &comm ):
    SubsetVariable( name ),
    d_iterator( iterator ),
    d_comm( comm )
{
}


AMP::Discretization::DOFManager::shared_ptr  MeshIteratorVariable::getSubsetDOF( AMP::Discretization::DOFManager::shared_ptr parentDOF )
{
    return parentDOF->subset( d_iterator, d_comm );
}


}
}

#endif

