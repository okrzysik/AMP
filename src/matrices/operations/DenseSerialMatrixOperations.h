#ifndef included_DenseSerialMatrixOperations_H_
#define included_DenseSerialMatrixOperations_H_

#include "AMP/matrices/operations/MatrixOperations.h"

namespace AMP::LinearAlgebra {

class DenseSerialMatrixOperations : public MatrixOperations
{

    /** \brief  Matrix-vector multiplication
     * \param[in]  x  The vector to multiply
     * \param[in]  A  The input matrix A
     * \param[out] y  The resulting vectory
     * \details  Compute \f$\mathbf{Ax} = \mathbf{y}\f$.
     */
    void mult( std::shared_ptr<const Vector> x,
               MatrixData const &A,
               std::shared_ptr<Vector> y ) override;

    /** \brief  Matrix transpose-vector multiplication
     * \param[in]  in  The vector to multiply
     * \param[in]  A  The input matrix A
     * \param[out] out The resulting vectory
     * \details  Compute \f$\mathbf{A}^T\mathbf{in} = \mathbf{out}\f$.
     */
    void multTranspose( std::shared_ptr<const Vector> in,
                        MatrixData const &A,
                        std::shared_ptr<Vector> out ) override;

    /** \brief  Scale the matrix by a scalar
     * \param[in] alpha  The value to scale by
     * \param[in,out] A The matrix A
     * \details  Compute \f$\mathbf{A} = \alpha\mathbf{A}\f$
     */
    void scale( AMP::Scalar alpha, MatrixData &A ) override;

    /** \brief  Scale the matrix by a scalar and diagonal matrix
     * \param[in] alpha  The value to scale by
     * \param[in] D  A vector representing the diagonal matrix
     * \param[in] A The input matrix A
     * \details  Compute \f$\mathbf{A} = \alpha\mathbf{D}\mathbf{A}\f$
     */
    void scale( AMP::Scalar, std::shared_ptr<const Vector>, MatrixData & ) override
    {
        AMP_ERROR( "Not implemented" );
    }

    /** \brief  Scale the matrix by a scalar and inverse of diagonal matrix
     * \param[in] alpha  The value to scale by
     * \param[in] D  A vector representing the diagonal matrix
     * \param[in] A The input matrix A
     * \details  Compute \f$\mathbf{A} = \alpha\mathbf{D}^{-1}\mathbf{A}\f$
     */
    void scaleInv( AMP::Scalar, std::shared_ptr<const Vector>, MatrixData & ) override
    {
        AMP_ERROR( "Not implemented" );
    }

    /** \brief  Compute the product of two matrices
     * \param[in] A  A multiplicand
     * \param[in] B  A multiplicand
     * \param[in] C  The product \f$\mathbf{AB}\f$.
     */
    void matMatMult( std::shared_ptr<MatrixData> A,
                     std::shared_ptr<MatrixData> B,
                     std::shared_ptr<MatrixData> C ) override;

    /** \brief  Compute the linear combination of two matrices
     * \param[in] alpha  scalar
     * \param[in] X matrix
     * \param[out] Y The output matrix
     * \details  Compute \f$\mathbf{THIS} = \alpha\mathbf{X} + \mathbf{THIS}\f$
     */
    void axpy( AMP::Scalar alpha, const MatrixData &X, MatrixData &Y ) override;

    /** \brief  Set the non-zeros of the matrix to a scalar
     * \param[in]  alpha  The value to set the non-zeros to
     * \param[out] A The input matrix A
     */
    void setScalar( AMP::Scalar alpha, MatrixData &A ) override;

    /** \brief  Set the non-zeros of the matrix to zero
     * \details  May not deallocate space.
     * \param[in] A The input matrix A
     */
    void zero( MatrixData &A ) override;

    /** \brief  Set the diagonal to the values in a vector
     * \param[in] in The values to set the diagonal to
     * \param[out] A The matrix to set
     */
    void setDiagonal( std::shared_ptr<const Vector> in, MatrixData &A ) override;

    /** \brief Extract the diagonal values into a vector
     * \param[in] in The values to set the diagonal to
     * \param[in] A The matrix to set
     */
    void extractDiagonal( MatrixData const &A, std::shared_ptr<Vector> buf ) override;

    /** \brief Extract the row sums into a vector
     * \param[in] A The matrix to read from
     * \param[out] buf Buffer to write row sums into
     */
    void getRowSums( MatrixData const &, std::shared_ptr<Vector> ) override
    {
        AMP_ERROR( "Not implemented" );
    }

    /** \brief Extract the absolute row sums into a vector
     * \param[in] A The matrix to read from
     * \param[out] buf Buffer to write row sums into
     */
    void getRowSumsAbsolute( MatrixData const &, std::shared_ptr<Vector> ) override
    {
        AMP_ERROR( "Not implemented" );
    }

    /** \brief  Set the matrix to the identity matrix
     * \param[out] A The matrix to set
     */
    void setIdentity( MatrixData &A ) override;

    /** \brief Compute the maximum row sum
     * \return  The L-infinity norm of the matrix
     * \param[in] X Data for the input matrix
     */
    AMP::Scalar LinfNorm( const MatrixData &X ) const override;

    /** \brief  Set <i>this</i> matrix with the same non-zero and distributed structure
     * as x and copy the coefficients
     * \param[in] x matrix data to copy from
     * \param[in] y matrix data to copy to
     */
    void copy( const MatrixData &x, MatrixData &y ) override;
};

} // namespace AMP::LinearAlgebra

#endif
