#include "AMP/vectors/VectorSelector.h"
#include "AMP/vectors/CommVariable.h"
#include "AMP/vectors/MeshVariable.h"
#include "AMP/vectors/StridedVariable.h"
#include "AMP/vectors/SubsetVariable.h"


namespace AMP::LinearAlgebra {


/********************************************************
 * VectorSelector                                        *
 ********************************************************/
VectorSelector::~VectorSelector() = default;
bool VectorSelector::isSelected( Vector::const_shared_ptr ) const { return true; }
AMP_MPI VectorSelector::communicator( Vector::const_shared_ptr p ) const { return p->getComm(); }
Vector::shared_ptr VectorSelector::subset( Vector::shared_ptr p ) const { return p; }
Vector::const_shared_ptr VectorSelector::subset( Vector::const_shared_ptr p ) const { return p; }


/********************************************************
 * VS_ByVariableName                                     *
 ********************************************************/
VS_ByVariableName::VS_ByVariableName( const std::string &name ) : d_VecName( name ) {}
bool VS_ByVariableName::isSelected( Vector::const_shared_ptr v ) const
{
    return v->getVariable()->getName() == d_VecName;
}
Vector::shared_ptr VS_ByVariableName::subset( Vector::shared_ptr vec ) const
{
    auto var = vec->getVariable();
    if ( var ) {
        if ( var->getName() == d_VecName )
            return vec;
    }
    return Vector::shared_ptr();
}
Vector::const_shared_ptr VS_ByVariableName::subset( Vector::const_shared_ptr vec ) const
{
    auto var = vec->getVariable();
    if ( var ) {
        if ( var->getName() == d_VecName )
            return vec;
    }
    return Vector::const_shared_ptr();
}


/********************************************************
 * VS_Stride                                             *
 ********************************************************/
VS_Stride::VS_Stride( size_t a, size_t b ) : d_Offset( a ), d_Stride( b ) {}
bool VS_Stride::isSelected( Vector::const_shared_ptr ) const { return true; }
Vector::shared_ptr VS_Stride::subset( Vector::shared_ptr p ) const
{
    auto variable =
        std::make_shared<StridedVariable>( p->getVariable()->getName(), d_Offset, d_Stride );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}
Vector::const_shared_ptr VS_Stride::subset( Vector::const_shared_ptr p ) const
{
    auto variable =
        std::make_shared<StridedVariable>( p->getVariable()->getName(), d_Offset, d_Stride );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}


/********************************************************
 * VS_Comm                                               *
 ********************************************************/
VS_Comm::VS_Comm( const AMP_MPI &comm )
{
    AMP_ASSERT( !comm.isNull() );
    d_comm = comm;
}
bool VS_Comm::isSelected( Vector::const_shared_ptr ) const { return true; }
AMP_MPI VS_Comm::communicator( Vector::const_shared_ptr p ) const
{
    return AMP_MPI::intersect( d_comm, p->getComm() );
}
Vector::shared_ptr VS_Comm::subset( Vector::shared_ptr p ) const
{
    auto variable =
        std::make_shared<CommVariable>( p->getVariable()->getName(), communicator( p ) );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}
Vector::const_shared_ptr VS_Comm::subset( Vector::const_shared_ptr p ) const
{
    auto variable =
        std::make_shared<CommVariable>( p->getVariable()->getName(), communicator( p ) );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}


/********************************************************
 * VS_Mesh                                               *
 ********************************************************/
VS_Mesh::VS_Mesh( AMP::Mesh::Mesh::shared_ptr mesh, bool useMeshComm )
    : d_useMeshComm( useMeshComm ), d_mesh( mesh )
{
}
bool VS_Mesh::isSelected( Vector::const_shared_ptr ) const { return true; }
AMP_MPI VS_Mesh::communicator( Vector::const_shared_ptr p ) const
{
    if ( d_useMeshComm ) {
        if ( d_mesh )
            return AMP_MPI::intersect( p->getComm(), d_mesh->getComm() );
        return AMP_MPI( AMP_COMM_WORLD );
    }
    return p->getComm();
}
Vector::shared_ptr VS_Mesh::subset( Vector::shared_ptr p ) const
{
    if ( !d_mesh )
        return Vector::shared_ptr();
    auto variable =
        std::make_shared<MeshVariable>( p->getVariable()->getName(), d_mesh, d_useMeshComm );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}
Vector::const_shared_ptr VS_Mesh::subset( Vector::const_shared_ptr p ) const
{
    if ( !d_mesh )
        return Vector::shared_ptr();
    auto variable =
        std::make_shared<MeshVariable>( p->getVariable()->getName(), d_mesh, d_useMeshComm );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}


/********************************************************
 * VS_MeshIterator                                       *
 ********************************************************/
VS_MeshIterator::VS_MeshIterator( const AMP::Mesh::MeshIterator &iterator,
                                  const AMP::AMP_MPI &comm )
    : d_comm( comm ), d_iterator( iterator )
{
}
bool VS_MeshIterator::isSelected( Vector::const_shared_ptr ) const { return true; }
Vector::shared_ptr VS_MeshIterator::subset( Vector::shared_ptr p ) const
{
    auto variable =
        std::make_shared<MeshIteratorVariable>( p->getVariable()->getName(), d_iterator, d_comm );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}
Vector::const_shared_ptr VS_MeshIterator::subset( Vector::const_shared_ptr p ) const
{
    auto variable =
        std::make_shared<MeshIteratorVariable>( p->getVariable()->getName(), d_iterator, d_comm );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}


} // namespace AMP::LinearAlgebra
