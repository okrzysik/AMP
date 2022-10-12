#ifndef included_AMP_DenseSerialMatrix
#define included_AMP_DenseSerialMatrix

#include "AMP/matrices/Matrix.h"
#include "AMP/vectors/Vector.h"
#include <memory>

namespace AMP::LinearAlgebra {


/** \class DenseSerialMatrix
 * \brief  An concrete class for dealing with dense serial matrices
 * \details  This is a concrete class that stores a dense local matrix.
 *    This is not a distributed matrix and requires that the comm is AMP_COMM_SELF.
 */
class DenseSerialMatrix : public Matrix
{
public:
    /** \brief Constructor
     * \param[in] params  Description of the matrix
     */
    explicit DenseSerialMatrix( std::shared_ptr<MatrixParameters> params );


    /** \brief Constructor
     * \param[in] params  MatrixData object associated with matrix
     */
    explicit DenseSerialMatrix( std::shared_ptr<MatrixData> data );

    DenseSerialMatrix( const DenseSerialMatrix & ) = delete;

    DenseSerialMatrix &operator=( const DenseSerialMatrix & ) = delete;

    /** \brief Destructor
     */
    virtual ~DenseSerialMatrix();

    //! Return the type of the matrix
    virtual std::string type() const override { return "DenseSerialMatrix"; }

    /** \brief  Matrix-vector multiplication
     * \param[in]  in  The vector to multiply
     * \param[out] out The resulting vectory
     * \details  Compute \f$\mathbf{Ain} = \mathbf{out}\f$.
     */
    void mult( AMP::LinearAlgebra::Vector::const_shared_ptr in,
               AMP::LinearAlgebra::Vector::shared_ptr out ) override;

    /** \brief  Matrix transpose-vector multiplication
     * \param[in]  in  The vector to multiply
     * \param[out] out The resulting vectory
     * \details  Compute \f$\mathbf{A}^T\mathbf{in} = \mathbf{out}\f$.
     */
    void multTranspose( AMP::LinearAlgebra::Vector::const_shared_ptr in,
                        AMP::LinearAlgebra::Vector::shared_ptr out ) override;


    /** \brief  Return a new matrix that is the transpose of this one
     * \return  A copy of this matrix transposed.
     */
    shared_ptr transpose() const override;

    /** \brief  Return a matrix with the same non-zero and distributed structure
     * \return  The new matrix
     */
    shared_ptr cloneMatrix() const override;

    /** \brief  Scale the matrix by a scalar
     * \param[in] alpha  The value to scale by
     * \details  Compute \f$\mathbf{A} = \alpha\mathbf{A}\f$
     */
    void scale( double alpha ) override;


    /** \brief  Compute the linear combination of two matrices
     * \param[in] alpha  scalar
     * \param[in] X matrix
     * \details  Compute \f$\mathbf{THIS} = \alpha\mathbf{X} + \mathbf{THIS}\f$
     */
    void axpy( double alpha, const Matrix &X ) override;


    /** \brief  Set the non-zeros of the matrix to a scalar
     * \param[in]  alpha  The value to set the non-zeros to
     */
    void setScalar( double alpha ) override;


    /** \brief  Set the non-zeros of the matrix to zero
     * \details  May not deallocate space.
     */
    void zero() override;


    /** \brief  Set the diagonal to the values in a vector
     * \param[in] in The values to set the diagonal to
     */
    void setDiagonal( Vector::const_shared_ptr in ) override;


    /** \brief  Set the matrix to the identity matrix
     */
    void setIdentity() override;

    /** \brief  Extract the diagonal from a matrix
     * \param[in]  buf  An optional vector to use as a buffer
     * \return  A vector of the diagonal values
     */
    Vector::shared_ptr
    extractDiagonal( Vector::shared_ptr buf = Vector::shared_ptr() ) const override;

    /** \brief Get a right vector( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$\mathbf{x}\f$ is a right
     * vector )
     * \return  A newly created right vector
     */
    Vector::shared_ptr getRightVector() const override;

    /** \brief Get a left vector( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$\mathbf{y}\f$ is a left
     * vector )
     * \return  A newly created left vector
     */
    Vector::shared_ptr getLeftVector() const override;

    /** \brief Compute the maximum column sum
     * \return  The L1 norm of the matrix
     */
    double L1Norm() const override;

protected:
    //! Unimplemented constructor
    DenseSerialMatrix();

    /** \brief  Multiply two matrices and store in a third
     * \param[in]  other_op  The other matrix to multiply
     * \param[out] result  The matrix to store the result
     */
    void multiply( shared_ptr other_op, shared_ptr &result ) override;
};
} // namespace AMP::LinearAlgebra


#endif
