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
bool VectorSelector::isSelected( const Vector & ) const { return true; }
AMP_MPI VectorSelector::communicator( const Vector &p ) const { return p.getComm(); }
std::shared_ptr<Vector> VectorSelector::subset( std::shared_ptr<Vector> p ) const { return p; }
std::shared_ptr<const Vector> VectorSelector::subset( std::shared_ptr<const Vector> p ) const
{
    return p;
}


/********************************************************
 * VS_ByVariableName                                     *
 ********************************************************/
VS_ByVariableName::VS_ByVariableName( const std::string &name ) : d_VecName( name ) {}
bool VS_ByVariableName::isSelected( const Vector &v ) const
{
    return v.getVariable()->getName() == d_VecName;
}
std::shared_ptr<Vector> VS_ByVariableName::subset( std::shared_ptr<Vector> vec ) const
{
    auto var = vec->getVariable();
    if ( var ) {
        if ( var->getName() == d_VecName )
            return vec;
    }
    return std::shared_ptr<Vector>();
}
std::shared_ptr<const Vector> VS_ByVariableName::subset( std::shared_ptr<const Vector> vec ) const
{
    auto var = vec->getVariable();
    if ( var ) {
        if ( var->getName() == d_VecName )
            return vec;
    }
    return std::shared_ptr<const Vector>();
}


/********************************************************
 * VS_Stride                                             *
 ********************************************************/
VS_Stride::VS_Stride( size_t a, size_t b ) : d_Offset( a ), d_Stride( b ) {}
bool VS_Stride::isSelected( const Vector & ) const { return true; }
std::shared_ptr<Vector> VS_Stride::subset( std::shared_ptr<Vector> p ) const
{
    auto variable =
        std::make_shared<StridedVariable>( p->getVariable()->getName(), d_Offset, d_Stride );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}
std::shared_ptr<const Vector> VS_Stride::subset( std::shared_ptr<const Vector> p ) const
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
bool VS_Comm::isSelected( const Vector & ) const { return true; }
AMP_MPI VS_Comm::communicator( const Vector &p ) const
{
    return AMP_MPI::intersect( d_comm, p.getComm() );
}
std::shared_ptr<Vector> VS_Comm::subset( std::shared_ptr<Vector> p ) const
{
    auto variable =
        std::make_shared<CommVariable>( p->getVariable()->getName(), communicator( *p ) );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}
std::shared_ptr<const Vector> VS_Comm::subset( std::shared_ptr<const Vector> p ) const
{
    auto variable =
        std::make_shared<CommVariable>( p->getVariable()->getName(), communicator( *p ) );
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
bool VS_Mesh::isSelected( const Vector & ) const { return true; }
AMP_MPI VS_Mesh::communicator( const Vector &p ) const
{
    if ( d_useMeshComm ) {
        if ( d_mesh )
            return AMP_MPI::intersect( p.getComm(), d_mesh->getComm() );
        return AMP_MPI( AMP_COMM_WORLD );
    }
    return p.getComm();
}
std::shared_ptr<Vector> VS_Mesh::subset( std::shared_ptr<Vector> p ) const
{
    if ( !d_mesh )
        return std::shared_ptr<Vector>();
    auto variable =
        std::make_shared<MeshVariable>( p->getVariable()->getName(), d_mesh, d_useMeshComm );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}
std::shared_ptr<const Vector> VS_Mesh::subset( std::shared_ptr<const Vector> p ) const
{
    if ( !d_mesh )
        return std::shared_ptr<Vector>();
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
bool VS_MeshIterator::isSelected( const Vector & ) const { return true; }
std::shared_ptr<Vector> VS_MeshIterator::subset( std::shared_ptr<Vector> p ) const
{
    auto variable =
        std::make_shared<MeshIteratorVariable>( p->getVariable()->getName(), d_iterator, d_comm );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}
std::shared_ptr<const Vector> VS_MeshIterator::subset( std::shared_ptr<const Vector> p ) const
{
    auto variable =
        std::make_shared<MeshIteratorVariable>( p->getVariable()->getName(), d_iterator, d_comm );
    auto vector = SubsetVariable::view( p, variable );
    return vector;
}


} // namespace AMP::LinearAlgebra
