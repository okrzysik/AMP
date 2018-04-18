#ifdef USE_AMP_VECTORS
#ifndef included_AMP_MatrixBuider
#define included_AMP_MatrixBuider

#include "AMP/matrices/Matrix.h"
#include "AMP/vectors/Vector.h"

#include <functional>
#include <string>


namespace AMP {
namespace LinearAlgebra {


/**
 * \brief  This function will create a matrix from two vectors
 * \details  This function is responsible for creating matrices given a left and a right vector
 * \param right     Vector that will be used to create the matrix
 *                  The right is x in the expression y = A*x.
 * \param left      Vector that will be used to create the matrix.
 *                  The left is y in the expression y = A*x.
 * \param type      Type of matrix to build:
 *                      auto: Automatically determined based on build (default)
 *                      ManagedPetscMatrix
 *                      ManagedEpetraMatrix
 *                      DenseSerialMatrix
 * \param getColumnIDs Function to provide the column indices given the row index.
 *                      If not provided, with will default to calling the getRowDOFs function on the
 *                      DOFManager associated with the left vector.
 */
AMP::LinearAlgebra::Matrix::shared_ptr
createMatrix( AMP::LinearAlgebra::Vector::shared_ptr right,
              AMP::LinearAlgebra::Vector::shared_ptr left,
              const std::string &type = "auto",
              std::function<std::vector<size_t>( size_t row )> getColumnIDs =
                  std::function<std::vector<size_t>( size_t )>() );

#if 0
/**
 * \brief  This function will create a matrix from two DOFManagers
 * \details  This function is responsible for creating matrices given left and right DOFManagers
 * \param right     Right DOFManager that determines the distribution of the columns
 * \param left      Left DOFManager that determines the distribution of the rows
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
    AMP::Discretization::DOFManager::shared_ptr right,
    AMP::Discretization::DOFManager::shared_ptr left,
    const std::string& type = "auto", 
    std::function<std::vector<size_t>(size_t)> getRow = std::function<std::vector<size_t>(size_t)>() );

#endif

} // namespace LinearAlgebra
} // namespace AMP

#endif
#endif
