#ifdef USE_AMP_VECTORS
#ifndef included_AMP_MatrixBuider
#define included_AMP_MatrixBuider

#include "matrices/Matrix.h"
#include "vectors/Vector.h"

#include <string>
#include <functional>


namespace AMP {
namespace LinearAlgebra {


/**
 * \brief  This function will create a matrix from two vectors
 * \details  This function is responsible for creating vectors from a DOFManager and variable.
 * \param right     Vector that will be used to create the matrix  The right is x in the expression y = A*x.
 * \param left      Vector that will be used to create the matrix.  The left is y in the expression y = A*x.
 * \param type      Type of matrix to build:
 *                      auto: Automatically determined based on build (default)
 *                      ManagedPetscMatrix
 *                      ManagedEpetraMatrix
 *                      DenseSerialMatrix
 * \param getRow    Function to provide the column indices given the row index.
 *                      If not provided, with will default to calling the getRowDOFs function on the
 *                      DOFManager associated with the left vector.
 */
AMP::LinearAlgebra::Matrix::shared_ptr createMatrix( 
    AMP::LinearAlgebra::Vector::shared_ptr right,
    AMP::LinearAlgebra::Vector::shared_ptr left,
    const std::string& type = "auto", 
    std::function<std::vector<size_t>(size_t)> getRow = std::function<std::vector<size_t>(size_t)>() );


} // LinearAlgebra namespace
} // AMP namespace

#endif
#endif
