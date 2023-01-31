#include "AMP/vectors/CommVariable.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/vectors/VectorSelector.h"

namespace AMP::LinearAlgebra {


CommVariable::CommVariable( const std::string &name, const AMP_MPI &comm ) : SubsetVariable( name )
{
    AMP_ASSERT( !comm.isNull() );
    d_comm = comm;
}


std::shared_ptr<AMP::Discretization::DOFManager>
CommVariable::getSubsetDOF( std::shared_ptr<AMP::Discretization::DOFManager> parentDOF ) const
{
    return parentDOF->subset( d_comm );
}


std::shared_ptr<VectorSelector> CommVariable::createVectorSelector() const
{
    return std::make_shared<VS_Comm>( d_comm );
}

} // namespace AMP::LinearAlgebra
