#include "vectors/VectorSelector.h"
#include "vectors/VectorSelector.h"
#include "vectors/SubsetVector.h"
#include "vectors/StridedVariable.h"
#include "vectors/MeshVariable.h"

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
* VS_Mesh                                               *
********************************************************/
#ifdef USE_AMP_MESH
VS_Mesh::VS_Mesh ( const std::string &name, AMP::Mesh::Mesh::shared_ptr mesh )
{
    d_Name = name;
    d_mesh = mesh;
}
Vector::shared_ptr  VS_Mesh::subset ( Vector::shared_ptr p ) const
{ 
    Variable::shared_ptr  variable ( new MeshVariable( d_Name, d_mesh ) );
    Vector::shared_ptr  vector = SubsetVector::view ( p, variable ); 
    return vector;
}
#endif


/********************************************************
* VS_MeshIterator                                       *
********************************************************/
#ifdef USE_AMP_MESH
VS_MeshIterator::VS_MeshIterator ( const std::string &name, const AMP::Mesh::MeshIterator &iterator ):
    d_iterator( iterator )
{
    d_Name = name;
}
Vector::shared_ptr  VS_MeshIterator::subset ( Vector::shared_ptr p ) const
{ 
    Variable::shared_ptr  variable ( new MeshIteratorVariable( d_Name, d_iterator ) );
    Vector::shared_ptr  vector = SubsetVector::view ( p, variable ); 
    return vector;
}
#endif


}
}

