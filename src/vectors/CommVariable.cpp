#include "CommVariable.h"
#include "AMP/discretization/DOF_Manager.h"

namespace AMP {
namespace LinearAlgebra {


CommVariable::CommVariable( const std::string &name, const AMP_MPI &comm ) : SubsetVariable( name )
{
    AMP_ASSERT( !comm.isNull() );
    d_comm = comm;
}


AMP::Discretization::DOFManager::shared_ptr
CommVariable::getSubsetDOF( std::shared_ptr<AMP::Discretization::DOFManager> parentDOF ) const
{
    return parentDOF->subset( d_comm );
}
} // namespace LinearAlgebra
} // namespace AMP
