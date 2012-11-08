#include "vectors/trilinos/NativeThyraVector.h"


namespace AMP {
namespace LinearAlgebra {

NativeThyraVector::NativeThyraVector ( VectorParameters::shared_ptr in_params ): 
    NativeVector(), 
    ThyraVector(), 
    VectorEngine()
{ 
    AMP_ERROR( "not implemented" );
}


NativeThyraVector::~NativeThyraVector ()
{
}


Vector::shared_ptr NativeThyraVector::cloneVector(const Variable::shared_ptr var ) const 
{ 
    AMP_ERROR( "not implemented" );
    return Vector::shared_ptr();
}

void NativeThyraVector::putRawData ( double *in )
{
    AMP_ERROR( "not implemented" );
}


}
}

