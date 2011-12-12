#include "VectorSubsetter.h"
#include "StridedVariable.h"


namespace AMP {
namespace LinearAlgebra {

VectorSubsetter::~VectorSubsetter () 
{
}

Vector::shared_ptr  VectorSubsetter::subset ( Vector::shared_ptr p ) 
{ 
    return p; 
}


VectorStriderSubsetter::VectorStriderSubsetter ( const std::string &n , size_t a , size_t b )
   : d_Offset ( a ) , d_Stride ( b ) , d_Name ( n ) 
{
}


VectorStriderSubsetter::~VectorStriderSubsetter () 
{
}


Vector::shared_ptr  VectorStriderSubsetter::subset ( Vector::shared_ptr p )
{ 
    Variable::shared_ptr  variable ( new StridedVariable( d_Name , d_Offset , d_Stride ) );
    Vector::shared_ptr  vector = SubsetVector::view ( p, variable ); 
    return vector;
}


VectorRandomAccessSubsetter::VectorRandomAccessSubsetter ( const std::string &n )
    : d_Name ( n )
    , d_RAIndexer ( new RandomAccessIndexer () )
{
}


Vector::shared_ptr   VectorRandomAccessSubsetter::subset ( Vector::shared_ptr p )
{
    AMP_ERROR("Not converted yet");
    return Vector::shared_ptr();
    // d_RAIndexer->castTo<RandomAccessIndexer>().finalize();
    // Variable *pNewVar = new RandomSubsetVariable ( d_Name , d_RAIndexer );
    // Variable::shared_ptr spNewVar ( pNewVar );
    // return SubsetVectorUsingSuperDOFs::view ( p , spNewVar );
}


}
}

