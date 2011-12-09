#ifdef USE_AMP_VECTORS
#ifndef included_AMP_MatrixBuider
#define included_AMP_MatrixBuider

#include <string>
#include "vectors/Vector.h"
#include "matrices/Matrix.h"


namespace AMP {
namespace LinearAlgebra {


/**
 * \brief  This function will create a matrix from two vectors
 * \details  This function is responsible for creating vectors from a DOFManager and variable.
 * \param operand  Variable that will be used to create the matrix
 * \param result   Variable that will be used to create the matrix
 */
AMP::LinearAlgebra::Matrix::shared_ptr  createMatrix( AMP::LinearAlgebra::Vector::shared_ptr operand , AMP::LinearAlgebra::Vector::shared_ptr result );


}
}

#endif
#endif

