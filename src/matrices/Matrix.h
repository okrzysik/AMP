#ifndef included_AMP_Matrix
#define included_AMP_Matrix

#include "AMP/matrices/MatrixParameters.h"
#include "AMP/utils/ParameterBase.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Vector.h"


namespace AMP {
namespace LinearAlgebra {


/** \class Matrix
 * \brief  An abstract interface for using and manipulating matrices
 * \details  There are several different varieties of distributed
 * memory matrices.  While most operations between the varieties
 * can be abstracted away from the user, some cannot.  For this
 * reason, most of the time, this class will suffice as the
 * way to interact with a matrix.  Matrix creation may require
 * use of one of the derived classes.
 */
class Matrix
{
public:
    //! Convenience typedef
    typedef AMP::shared_ptr<Matrix> shared_ptr;
    typedef AMP::shared_ptr<const Matrix> const_shared_ptr;

    /** \brief Constructor
     * \param[in] params  Description of the matrix
     */
    explicit Matrix( MatrixParameters::shared_ptr params );

    /** \brief Destructor
     */
    virtual ~Matrix();

    /** \brief  Matrix-vector multiplication
     * \param[in]  in  The vector to multiply
     * \param[out] out The resulting vectory
     * \details  Compute \f$\mathbf{Ain} = \mathbf{out}\f$.
     */
    virtual void mult( AMP::LinearAlgebra::Vector::const_shared_ptr in,
                       AMP::LinearAlgebra::Vector::shared_ptr out ) = 0;

    /** \brief  Matrix transpose-vector multiplication
     * \param[in]  in  The vector to multiply
     * \param[out] out The resulting vectory
     * \details  Compute \f$\mathbf{A}^T\mathbf{in} = \mathbf{out}\f$.
     */
    virtual void multTranspose( AMP::LinearAlgebra::Vector::const_shared_ptr in,
                                AMP::LinearAlgebra::Vector::shared_ptr out ) = 0;


    /** \brief  Return a new matrix that is the transpose of this one
     * \return  A copy of this matrix transposed.
     */
    virtual shared_ptr transpose() const;

    /** \brief  Return a matrix with the same non-zero and distributed structure
     * \return  The new matrix
     */
    virtual shared_ptr cloneMatrix() const = 0;

    /** \brief  Scale the matrix by a scalar
     * \param[in] alpha  The value to scale by
     * \details  Compute \f$\mathbf{A} = \alpha\mathbf{A}\f$
     */
    virtual void scale( double alpha ) = 0;


    /** \brief  Compute the product of two matrices
     * \param[in] A  A multiplicand
     * \param[in] B  A multiplicand
     * \return The product \f$\mathbf{AB}\f$.
     */
    static shared_ptr matMultiply( shared_ptr A, shared_ptr B );

    /** \brief  Compute the linear combination of two matrices
     * \param[in] alpha  scalar
     * \param[in] X matrix
     * \details  Compute \f$\mathbf{THIS} = \alpha\mathbf{X} + \mathbf{THIS}\f$
     */
    virtual void axpy( double alpha, const Matrix &X ) = 0;

    /** \brief  Compute the linear combination of two matrices
     * \param[in] alpha  scalar
     * \param[in] X matrix
     * \details  Compute \f$\mathbf{THIS} = \alpha\mathbf{X} + \mathbf{THIS}\f$
     */
    void axpy( double alpha, Matrix::const_shared_ptr X );


    /** \brief  Add values to those in the matrix
     * \param[in] num_rows The number of rows represented in values
     * \param[in] num_cols The number of cols represented in values
     * \param[in] rows  The row ids of values
     * \param[in] cols  The column ids of values
     * \param[in] values  The values to add to the matrix (row-major ordering)
     * \details  This method may fail if the matrix has not
     * allocated a particular (row,col) specified, depending
     * on the actual subclass of matrix used.
     */
    virtual void addValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) = 0;

    /** \brief  Set values in the matrix
     * \param[in] num_rows The number of rows represented in values
     * \param[in] num_cols The number of cols represented in values
     * \param[in] rows  The row ids of values
     * \param[in] cols  The column ids of values
     * \param[in] values  The values to set to the matrix (row-major ordering)
     * \details  This method may fail if the matrix has not
     * allocated a particular (row,col) specified, depending
     * on the actual subclass of matrix used.
     */
    virtual void setValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) = 0;

    /** \brief  Get values in the matrix
     * \param[in] num_rows The number of rows represented in values
     * \param[in] num_cols The number of cols represented in values
     * \param[in] rows  The row ids of values
     * \param[in] cols  The column ids of values
     * \param[out] values  The values to get from the matrix (row-major ordering)
     * \details  This method will return zero for any entries that
     *   have not been allocated or are not ghosts on the current processor.
     */
    virtual void getValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) const = 0;


    /** \brief  Add values to those in the matrix
     * \param[in] row  The row id of value
     * \param[in] col  The column id of value
     * \param[in] value  The value to add to the matrix
     * \details  This method may fail if the matrix has not
     * allocated a particular (row,col) specified, depending
     * on the actual subclass of matrix used.
     */
    virtual void addValueByGlobalID( size_t row, size_t col, double value );

    /** \brief  Set values in the matrix
     * \param[in] row  The row id of value
     * \param[in] col  The column id of value
     * \param[in] value  The value to set to the matrix
     * \details  This method may fail if the matrix has not
     * allocated a particular (row,col) specified, depending
     * on the actual subclass of matrix used.
     */
    virtual void setValueByGlobalID( size_t row, size_t col, double value );

    /** \brief  Get values in the matrix
     * \param[in] row  The row id of value
     * \param[in] col  The column id of value
     * \details  This method will return zero for any values that have not been allocated.
     */
    virtual double getValueByGlobalID( size_t row, size_t col ) const;


    /** \brief  Set the non-zeros of the matrix to a scalar
     * \param[in]  alpha  The value to set the non-zeros to
     */
    virtual void setScalar( double alpha ) = 0;

    /** \brief  Set the non-zeros of the matrix to zero
     * \details  May not deallocate space.
     */
    virtual void zero() = 0;

    /** \brief  Retrieve a row of the matrix in compressed format
     * \param[in]  row Which row
     * \param[out] cols  The column ids of the returned values
     * \param[out] values  The values in the row
     */
    virtual void getRowByGlobalID( size_t row,
                                   std::vector<size_t> &cols,
                                   std::vector<double> &values ) const = 0;

    /** \brief  Given a row, retrieve the non-zero column indices of the matrix in compressed format
     * \param[in]  row Which row
     */
    virtual std::vector<size_t> getColumnIDs( size_t row ) const = 0;

    /** \brief  Set the diagonal to the values in a vector
     * \param[in] in The values to set the diagonal to
     */
    virtual void setDiagonal( Vector::const_shared_ptr in ) = 0;

    /** \brief  Set the matrix to the identity matrix
     */
    virtual void setIdentity() = 0;

    /** \brief  Perform communication to ensure values in the
     * matrix are the same across cores.
     */
    virtual void makeConsistent() = 0;

    /** \brief  Get the number of local rows in the matrix
     * \return  The number of local rows
     */
    virtual size_t numLocalRows() const;

    /** \brief  Get the number of global rows in the matrix
     * \return  The number of global rows
     */
    virtual size_t numGlobalRows() const;

    /** \brief  Get the number of local columns in the matrix
     * \return  The number of local columns
     */
    virtual size_t numLocalColumns() const;

    /** \brief  Get the number of global columns in the matrix
     * \return  The number of global columns
     */
    virtual size_t numGlobalColumns() const;

    /** \brief  Get the global id of the beginning row
     * \return  beginning global row id
     */
    virtual size_t beginRow() const;

    /** \brief  Get the global id of the ending row
     * \return  ending global row id
     */
    virtual size_t endRow() const;

    /** \brief  Extract the diagonal from a matrix
     * \param[in]  buf  An optional vector to use as a buffer
     * \return  A vector of the diagonal values
     */
    virtual Vector::shared_ptr
    extractDiagonal( Vector::shared_ptr buf = Vector::shared_ptr() ) const = 0;

    /** \brief Get a right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$\mathbf{x}\f$ is a right
     * vector )
     * \return  A newly created right vector
     */
    virtual Vector::shared_ptr getRightVector() const = 0;

    /** \brief Get a left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$\mathbf{y}\f$ is a left
     * vector )
     * \return  A newly created left vector
     */
    virtual Vector::shared_ptr getLeftVector() const = 0;

    /** \brief Get the DOFManager associated with a right vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$\mathbf{x}\f$
     * is a right vector )
     * \return  The DOFManager associated with a right vector
     */
    virtual Discretization::DOFManager::shared_ptr getRightDOFManager() const = 0;

    /** \brief Get the DOFManager associated with a left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$\mathbf{y}\f$ is
     * a left vector )
     * \return  The DOFManager associated with a left vector
     */
    virtual Discretization::DOFManager::shared_ptr getLeftDOFManager() const = 0;

    /** \brief Compute the maximum column sum
     * \return  The L1 norm of the matrix
     */
    virtual double L1Norm() const = 0;


protected:
    //! Protected constructor
    Matrix();

    //! Protected copy constructor
    Matrix( const Matrix & );

    //! Protected assignment operator
    Matrix &operator=( const Matrix & ) = delete;


protected:
    /** \brief  Multiply two matrices and store in a third
     *    result = this * other_op
     * \param[in]  other_op  The other matrix to multiply
     * \param[out] result  The matrix to store the result
     */
    virtual void multiply( shared_ptr other_op, shared_ptr &result ) = 0;

    //!  Communicator for the matrix
    AMP_MPI d_comm;
};


//! Stream operator
std::ostream &operator<<( std::ostream &out, const Matrix::shared_ptr );
//! Stream operator
std::ostream &operator<<( std::ostream &out, const Matrix & );
} // namespace LinearAlgebra
} // namespace AMP

#include "Matrix.inline.h"
#endif
