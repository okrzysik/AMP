#include "AMP/vectors/MeshVariable.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/vectors/VectorSelector.h"


namespace AMP::LinearAlgebra {


/****************************************************************
 * MeshVariable                                                  *
 ****************************************************************/
MeshVariable::MeshVariable( const std::string &name,
                            std::shared_ptr<AMP::Mesh::Mesh> mesh,
                            bool useMeshComm )
    : SubsetVariable( name )
{
    AMP_ASSERT( mesh );
    d_mesh        = mesh;
    d_useMeshComm = useMeshComm;
}
std::shared_ptr<AMP::Discretization::DOFManager>
MeshVariable::getSubsetDOF( std::shared_ptr<AMP::Discretization::DOFManager> parentDOF ) const
{
    return parentDOF->subset( d_mesh, d_useMeshComm );
}
std::shared_ptr<VectorSelector> MeshVariable::createVectorSelector() const
{
    return std::make_shared<VS_Mesh>( d_mesh, d_useMeshComm );
}


/****************************************************************
 * MeshVariable                                                  *
 ****************************************************************/
MeshIteratorVariable::MeshIteratorVariable( const std::string &name,
                                            const AMP::Mesh::MeshIterator &iterator,
                                            const AMP_MPI &comm )
    : SubsetVariable( name ), d_comm( std::move( comm ) ), d_iterator( iterator )
{
}
std::shared_ptr<AMP::Discretization::DOFManager> MeshIteratorVariable::getSubsetDOF(
    std::shared_ptr<AMP::Discretization::DOFManager> parentDOF ) const
{
    return parentDOF->subset( d_iterator, d_comm );
}
std::shared_ptr<VectorSelector> MeshIteratorVariable::createVectorSelector() const
{
    return std::make_shared<VS_MeshIterator>( d_iterator, d_comm );
}


} // namespace AMP::LinearAlgebra
