
#include "utils/ParameterBase.h"
#include "Matrix.h"

namespace AMP {
namespace LinearAlgebra {

  Matrix::shared_ptr  Matrix::matMultiply ( shared_ptr A , shared_ptr B )
  {
    shared_ptr retVal;
    A->multiply ( B , retVal );
    return retVal;
  }

}
}

