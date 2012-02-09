#include "StridedVariable.h"

namespace AMP {
namespace LinearAlgebra {


StridedVariable::StridedVariable ( const std::string &name, size_t offset, size_t stride ):
    SubsetVariable ( name )
{
    d_offset = offset;
    d_stride = stride;
}


AMP::Discretization::DOFManager::shared_ptr  StridedVariable::getSubsetDOF( AMP::Discretization::DOFManager::shared_ptr parentDOF )
{ 
    std::vector<size_t> dofs;
    dofs.reserve(parentDOF->numLocalDOF());
    size_t offset2 = d_stride - d_offset%d_stride;  // Change the sign of the offset
    size_t i = 0;
    for (size_t DOF=parentDOF->beginDOF(); DOF<parentDOF->endDOF(); DOF++) {
        if ( (DOF+offset2)%d_stride==0 ) {
            dofs.push_back(DOF);
            i++;
        }
    }
    AMP::Mesh::MeshIterator iterator = parentDOF->getIterator();
    AMP::Discretization::DOFManager::shared_ptr  subsetDOF;
    subsetDOF = AMP::Discretization::subsetDOFManager::create( parentDOF, dofs, iterator, parentDOF->getComm() );
    return subsetDOF;
}


}
}
