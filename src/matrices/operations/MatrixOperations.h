#ifndef included_MatrixOperations_H_
#define included_MatrixOperations_H_

#include "AMP/vectors/Scalar.h"

namespace AMP::LinearAlgebra {

class Vector;
class MatrixData;

class MatrixOperations
{
public:
    //! Destructor
    virtual ~MatrixOperations() {}

    /** \brief  Matrix-vector multiplication
     * \param[in]  in  The vector to multiply
     * \param[out] out The resulting vectory
     * \details  Compute \f$\mathbf{Ax} = \mathbf{y}\f$.
     */
    virtual void
    mult( std::shared_ptr<const Vector> x, MatrixData const &A, std::shared_ptr<Vector> y ) = 0;

    /** \brief  Matrix transpose-vector multiplication
     * \param[in]  in  The vector to multiply
     * \param[out] out The resulting vectory
     * \details  Compute \f$\mathbf{A}^T\mathbf{in} = \mathbf{out}\f$.
     */
    virtual void multTranspose( std::shared_ptr<const Vector> in,
                                MatrixData const &A,
                                std::shared_ptr<Vector> out ) = 0;

    /** \brief  Scale the matrix by a scalar
     * \param[in] alpha  The value to scale by
     * \details  Compute \f$\mathbf{A} = \alpha\mathbf{A}\f$
     */
    virtual void scale( AMP::Scalar alpha, MatrixData &A ) = 0;

    /** \brief  Compute the product of two matrices
     * \param[in] A  A multiplicand
     * \param[in] B  A multiplicand
     * \return The product \f$\mathbf{AB}\f$.
     */
    virtual void matMultiply( MatrixData const &A, MatrixData const &B, MatrixData &C ) = 0;

    /** \brief  Compute the linear combination of two matrices
     * \param[in] alpha  scalar
     * \param[in] X matrix
     * \details  Compute \f$\mathbf{THIS} = \alpha\mathbf{X} + \mathbf{THIS}\f$
     */
    virtual void axpy( AMP::Scalar alpha, const MatrixData &X, MatrixData &Y ) = 0;

    /** \brief  Set the non-zeros of the matrix to a scalar
     * \param[in]  alpha  The value to set the non-zeros to
     */
    virtual void setScalar( AMP::Scalar alpha, MatrixData &A ) = 0;

    /** \brief  Set the non-zeros of the matrix to zero
     * \details  May not deallocate space.
     */
    virtual void zero( MatrixData &A ) = 0;

    /** \brief  Set the diagonal to the values in a vector
     * \param[in] in The values to set the diagonal to
     */
    virtual void setDiagonal( std::shared_ptr<const Vector> in, MatrixData &A ) = 0;

    /** \brief  Set the matrix to the identity matrix
     */
    virtual void setIdentity( MatrixData &A ) = 0;

    /** \brief Compute the maximum column sum
     * \return  The L1 norm of the matrix
     */
    virtual AMP::Scalar L1Norm( const MatrixData &X ) const = 0;
};

} // namespace AMP::LinearAlgebra

#endif // included_MatrixOperations_H_
