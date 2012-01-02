#ifdef USE_AMP_MESH
#include "MeshVariable.h"

namespace AMP {
namespace LinearAlgebra {


MeshVariable::MeshVariable ( const std::string &name, AMP::Mesh::Mesh::shared_ptr mesh ):
    SubsetVariable ( name )
{
    d_mesh = mesh;
}


boost::shared_ptr<AMP::Discretization::subsetDOFManager>  MeshVariable::getSubsetDOF( AMP::Discretization::DOFManager::shared_ptr parentDOF )
{
    return boost::static_pointer_cast<AMP::Discretization::subsetDOFManager>( parentDOF->subset( d_mesh ) );
}


MeshIteratorVariable::MeshIteratorVariable ( const std::string &name, const AMP::Mesh::MeshIterator &iterator ):
    SubsetVariable( name ),
    d_iterator( iterator )
{
}


boost::shared_ptr<AMP::Discretization::subsetDOFManager>  MeshIteratorVariable::getSubsetDOF( AMP::Discretization::DOFManager::shared_ptr parentDOF )
{
    return boost::static_pointer_cast<AMP::Discretization::subsetDOFManager>( parentDOF->subset( d_iterator ) );
}


}
}

#endif

