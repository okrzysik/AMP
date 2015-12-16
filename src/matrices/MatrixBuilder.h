#ifdef USE_AMP_VECTORS
#ifndef included_AMP_MatrixBuider
#define included_AMP_MatrixBuider

#include "matrices/Matrix.h"
#include "vectors/Vector.h"
#include <string>


namespace AMP {
namespace LinearAlgebra {


/**
 * \brief  This function will create a matrix from two vectors
 * \details  This function is responsible for creating vectors from a DOFManager and variable.
 * \param operand  Vector that will be used to create the matrix.  The operand vector is the right
 *                 vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector ) and
 *                 determines the number of rows.
 * \param result   Vector that will be used to create the matrix  The result vector is the left
 *                 vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector ) and
 *                 determines the number of columns.
 * \param type     Type of matrix to build:
 *                 0: Automatically determined based on build (default)
 *                 1: ManagedPetscMatrix
 *                 2: DenseSerialMatrix
 */
AMP::LinearAlgebra::Matrix::shared_ptr createMatrix( AMP::LinearAlgebra::Vector::shared_ptr operand,
                                                     AMP::LinearAlgebra::Vector::shared_ptr result,
                                                     int type = 0 );
}
}

#endif
#endif
