#include "AMP/vectors/VectorSelector.h"

#include "AMP/vectors/CommVariable.h"
#include "AMP/vectors/MeshVariable.h"
#include "AMP/vectors/StridedVariable.h"
#include "AMP/vectors/SubsetVector.h"
#include <utility>


namespace AMP {
namespace LinearAlgebra {


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


/********************************************************
 * VS_Stride                                             *
 ********************************************************/
VS_Stride::VS_Stride( size_t a, size_t b ) : d_Offset( a ), d_Stride( b ) {}
Vector::shared_ptr VS_Stride::subset( Vector::shared_ptr p ) const
{
    auto variable =
        AMP::make_shared<StridedVariable>( p->getVariable()->getName(), d_Offset, d_Stride );
    auto vector = SubsetVector::view( p, variable );
    return vector;
}
Vector::const_shared_ptr VS_Stride::subset( Vector::const_shared_ptr p ) const
{
    auto variable =
        AMP::make_shared<StridedVariable>( p->getVariable()->getName(), d_Offset, d_Stride );
    auto vector = SubsetVector::view( p, variable );
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
AMP_MPI VS_Comm::communicator( Vector::const_shared_ptr p ) const
{
    return AMP_MPI::intersect( d_comm, p->getComm() );
}
Vector::shared_ptr VS_Comm::subset( Vector::shared_ptr p ) const
{
    auto variable =
        AMP::make_shared<CommVariable>( p->getVariable()->getName(), communicator( p ) );
    auto vector = SubsetVector::view( p, variable );
    return vector;
}
Vector::const_shared_ptr VS_Comm::subset( Vector::const_shared_ptr p ) const
{
    auto variable =
        AMP::make_shared<CommVariable>( p->getVariable()->getName(), communicator( p ) );
    auto vector = SubsetVector::view( p, variable );
    return vector;
}


/********************************************************
 * VS_Mesh                                               *
 ********************************************************/
#ifdef USE_AMP_MESH
VS_Mesh::VS_Mesh( AMP::Mesh::Mesh::shared_ptr mesh, bool useMeshComm )
{
    d_mesh        = mesh;
    d_useMeshComm = useMeshComm;
}
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
    if ( d_mesh == nullptr )
        return Vector::shared_ptr();
    auto variable =
        AMP::make_shared<MeshVariable>( p->getVariable()->getName(), d_mesh, d_useMeshComm );
    auto vector = SubsetVector::view( p, variable );
    return vector;
}
Vector::const_shared_ptr VS_Mesh::subset( Vector::const_shared_ptr p ) const
{
    if ( d_mesh == nullptr )
        return Vector::shared_ptr();
    auto variable =
        AMP::make_shared<MeshVariable>( p->getVariable()->getName(), d_mesh, d_useMeshComm );
    auto vector = SubsetVector::view( p, variable );
    return vector;
}
#endif


/********************************************************
 * VS_MeshIterator                                       *
 ********************************************************/
#ifdef USE_AMP_MESH
VS_MeshIterator::VS_MeshIterator( const AMP::Mesh::MeshIterator &iterator,
                                  const AMP::AMP_MPI &comm )
    : d_comm( comm ), d_iterator( iterator )
{
}
Vector::shared_ptr VS_MeshIterator::subset( Vector::shared_ptr p ) const
{
    auto variable =
        AMP::make_shared<MeshIteratorVariable>( p->getVariable()->getName(), d_iterator, d_comm );
    auto vector = SubsetVector::view( p, variable );
    return vector;
}
Vector::const_shared_ptr VS_MeshIterator::subset( Vector::const_shared_ptr p ) const
{
    auto variable =
        AMP::make_shared<MeshIteratorVariable>( p->getVariable()->getName(), d_iterator, d_comm );
    auto vector = SubsetVector::view( p, variable );
    return vector;
}
#endif
} // namespace LinearAlgebra
} // namespace AMP
