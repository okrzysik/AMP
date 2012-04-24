#include "vectors/VectorSelector.h"
#include "vectors/VectorSelector.h"
#include "vectors/SubsetVector.h"
#include "vectors/StridedVariable.h"
#include "vectors/MeshVariable.h"
#include "vectors/CommVariable.h"

namespace AMP {
namespace LinearAlgebra {


/********************************************************
* VectorSelector                                        *
********************************************************/
VectorSelector::~VectorSelector ()
{
}
bool  VectorSelector::isSelected ( Vector::const_shared_ptr ) const
{
    return true;
}
AMP_MPI  VectorSelector::communicator ( Vector::const_shared_ptr p ) const
{
    return p->getComm();
}
Vector::shared_ptr  VectorSelector::subset ( Vector::shared_ptr p ) const
{
    return p;
}

  

/********************************************************
* VS_ByVariableName                                     *
********************************************************/
VS_ByVariableName::VS_ByVariableName ( std::string  n )
    : d_VecName ( n ) 
{
}
bool   VS_ByVariableName::isSelected ( Vector::const_shared_ptr v ) const
{
    return v->getVariable()->getName() == d_VecName;
}


/********************************************************
* VS_Stride                                             *
********************************************************/
VS_Stride::VS_Stride ( const std::string &n , size_t a , size_t b ) : 
    d_Offset ( a ),
    d_Stride ( b ),
    d_Name( n ) 
{
}
Vector::shared_ptr  VS_Stride::subset ( Vector::shared_ptr p ) const
{ 
    Variable::shared_ptr  variable ( new StridedVariable( d_Name , d_Offset , d_Stride ) );
    Vector::shared_ptr  vector = SubsetVector::view ( p, variable ); 
    return vector;
}


/********************************************************
* VS_Comm                                               *
********************************************************/
VS_Comm::VS_Comm ( const std::string &name, AMP_MPI comm )
{
    d_Name = name;
    AMP_ASSERT(!comm.isNull());
    d_comm = comm;
}
AMP_MPI  VS_Comm::communicator ( Vector::const_shared_ptr p ) const
{
    return AMP_MPI::intersect( d_comm, p->getComm() );
}
Vector::shared_ptr  VS_Comm::subset ( Vector::shared_ptr p ) const
{ 
    Variable::shared_ptr  variable ( new CommVariable( d_Name, communicator(p) ) );
    Vector::shared_ptr  vector = SubsetVector::view ( p, variable ); 
    return vector;
}


/********************************************************
* VS_Mesh                                               *
********************************************************/
#ifdef USE_AMP_MESH
VS_Mesh::VS_Mesh ( const std::string &name, AMP::Mesh::Mesh::shared_ptr mesh, bool useMeshComm )
{
    AMP_ASSERT(mesh.get()!=NULL);
    d_Name = name;
    d_mesh = mesh;
    d_useMeshComm = useMeshComm;
}
AMP_MPI  VS_Mesh::communicator ( Vector::const_shared_ptr p ) const
{
    if ( d_useMeshComm )
        return AMP_MPI::intersect( p->getComm(), d_mesh->getComm() );
    return p->getComm();
}
Vector::shared_ptr  VS_Mesh::subset ( Vector::shared_ptr p ) const
{ 
    Variable::shared_ptr  variable ( new MeshVariable( d_Name, d_mesh, d_useMeshComm ) );
    Vector::shared_ptr  vector = SubsetVector::view ( p, variable ); 
    return vector;
}
#endif


/********************************************************
* VS_MeshIterator                                       *
********************************************************/
#ifdef USE_AMP_MESH
VS_MeshIterator::VS_MeshIterator ( const std::string &name, const AMP::Mesh::MeshIterator &iterator, const AMP::AMP_MPI &comm ):
    d_iterator( iterator ),
    d_comm( comm )
{
    d_Name = name;
}
Vector::shared_ptr  VS_MeshIterator::subset ( Vector::shared_ptr p ) const
{ 
    Variable::shared_ptr  variable ( new MeshIteratorVariable( d_Name, d_iterator, d_comm ) );
    Vector::shared_ptr  vector = SubsetVector::view ( p, variable ); 
    return vector;
}
#endif


}
}

