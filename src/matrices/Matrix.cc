
#include "utils/ParameterBase.h"
#include "Matrix.h"

namespace AMP {
namespace LinearAlgebra {

Matrix::shared_ptr  Matrix::matMultiply ( shared_ptr A , shared_ptr B )
{
    if ( A->numColumns() != B->numRows() )
        AMP_ERROR( "Inner matrix dimensions must agree" );
    shared_ptr retVal;
    A->multiply ( B , retVal );
    return retVal;
}

}
}

