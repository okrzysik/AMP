#include "CommVariable.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
namespace LinearAlgebra {


CommVariable::CommVariable ( const std::string &name, AMP_MPI comm ):
    SubsetVariable ( name )
{
    AMP_ASSERT(!comm.isNull());
    d_comm = comm;
}


AMP::Discretization::DOFManager::shared_ptr  CommVariable::getSubsetDOF( AMP::Discretization::DOFManager::shared_ptr parentDOF ) const
{
    return parentDOF->subset( d_comm );
}


}
}


