#ifndef included_AMP_DenseSerialMatrix
#define included_AMP_DenseSerialMatrix

#include "matrices/Matrix.h"
#include "utils/shared_ptr.h"
#include "vectors/Vector.h"

namespace AMP {
namespace LinearAlgebra {


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
    explicit DenseSerialMatrix( MatrixParameters::shared_ptr params );

    /** \brief Destructor
      */
    virtual ~DenseSerialMatrix();

    /** \brief  Matrix-vector multiplication
      * \param[in]  in  The vector to multiply
      * \param[out] out The resulting vectory
      * \details  Compute \f$\mathbf{Ain} = \mathbf{out}\f$.
      */
    virtual void mult( AMP::LinearAlgebra::Vector::const_shared_ptr in,
                       AMP::LinearAlgebra::Vector::shared_ptr out );

    /** \brief  Matrix transpose-vector multiplication
      * \param[in]  in  The vector to multiply
      * \param[out] out The resulting vectory
      * \details  Compute \f$\mathbf{A}^T\mathbf{in} = \mathbf{out}\f$.
      */
    virtual void multTranspose( AMP::LinearAlgebra::Vector::const_shared_ptr in,
                                AMP::LinearAlgebra::Vector::shared_ptr out );


    /** \brief  Return a new matrix that is the transpose of this one
      * \return  A copy of this matrix transposed.
      */
    virtual shared_ptr transpose() const;

    /** \brief  Return a matrix with the same non-zero and distributed structure
      * \return  The new matrix
      */
    virtual shared_ptr cloneMatrix() const;

    /** \brief  Scale the matrix by a scalar
      * \param[in] alpha  The value to scale by
      * \details  Compute \f$\mathbf{A} = \alpha\mathbf{A}\f$
      */
    virtual void scale( double alpha );


    /** \brief  Compute the linear combination of two matrices
      * \param[in] alpha  scalar
      * \param[in] X matrix
      * \details  Compute \f$\mathbf{THIS} = \alpha\mathbf{X} + \mathbf{THIS}\f$
      */
    virtual void axpy( double alpha, const Matrix &X );


    /** \brief  Add values to those in the matrix
      * \param[in] num_rows The number of rows represented in values
      * \param[in] num_cols The number of cols represented in values
      * \param[in] rows  The row ids of values
      * \param[in] cols  The column ids of values
      * \param[in] values  The values to add to the matrix
      * \details  This method may fail if the matrix has not
      * allocated a particular(row,col) specified, depending
      * on the actual subclass of matrix used.
      */
    virtual void
    addValuesByGlobalID( size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values );

    /** \brief  Set values in the matrix
      * \param[in] num_rows The number of rows represented in values
      * \param[in] num_cols The number of cols represented in values
      * \param[in] rows  The row ids of values
      * \param[in] cols  The column ids of values
      * \param[in] values  The values to set to the matrix
      * \details  This method may fail if the matrix has not
      * allocated a particular(row,col) specified, depending
      * on the actual subclass of matrix used.
      */
    virtual void
    setValuesByGlobalID( size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values );

    /** \brief  Get values in the matrix
      * \param[in] num_rows The number of rows represented in values
      * \param[in] num_cols The number of cols represented in values
      * \param[in] rows  The row ids of values
      * \param[in] cols  The column ids of values
      * \param[in] values  The values to get from the matrix (row-major ordering)
      * \details  This method will return zero for any entries that
      *   have not been allocated or are not ghosts on the current processor.
      */
    virtual void
    getValuesByGlobalID( size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) const;


    /** \brief  Add values to those in the matrix
      * \param[in] row  The row id of value
      * \param[in] col  The column id of value
      * \param[in] value  The value to add to the matrix
      * \details  This method may fail if the matrix has not
      * allocated a particular(row,col) specified, depending
      * on the actual subclass of matrix used.
      */
    virtual void addValueByGlobalID( size_t row, size_t col, double value );

    /** \brief  Set values in the matrix
      * \param[in] row  The row id of value
      * \param[in] col  The column id of value
      * \param[in] value  The value to set to the matrix
      * \details  This method may fail if the matrix has not
      * allocated a particular(row,col) specified, depending
      * on the actual subclass of matrix used.
      */
    virtual void setValueByGlobalID( size_t row, size_t col, double value );

    /** \brief  Set values in the matrix
      * \param[in] row  The row id of value
      * \param[in] col  The column id of value
      * \details  This method may fail if the matrix has not
      * allocated a particular(row,col) specified, depending
      * on the actual subclass of matrix used.
      */
    virtual double getValueByGlobalID( size_t row, size_t col ) const;


    /** \brief  Set the non-zeros of the matrix to a scalar
      * \param[in]  alpha  The value to set the non-zeros to
      */
    virtual void setScalar( double alpha );


    /** \brief  Retrieve a row of the matrix in compressed format
      * \param[in]  row Which row
      * \param[out] cols  The column ids of the returned values
      * \param[out] values  The values in the row
      */
    virtual void
    getRowByGlobalID( size_t row, std::vector<size_t> &cols, std::vector<double> &values ) const;

    /** \brief  Set the diagonal to the values in a vector
      * \param[in] in The values to set the diagonal to
      */
    virtual void setDiagonal( Vector::const_shared_ptr in );

    /** \brief  Set the matrix to the identity matrix
     */
    virtual void setIdentity();

    /** \brief  Set the non-zeros of the matrix to zero
      * \details  May not deallocate space.
      */
    virtual void zero();

    /** \brief  Perform communication to ensure values in the
      * matrix are the same across cores.
      */
    virtual void makeConsistent() {}

    /** \brief  Extract the diagonal from a matrix
      * \param[in]  buf  An optional vector to use as a buffer
      * \return  A vector of the diagonal values
      */
    virtual Vector::shared_ptr
    extractDiagonal( Vector::shared_ptr buf = Vector::shared_ptr() ) const;

    /** \brief Get a right vector( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$\mathbf{x}\f$ is a right
     * vector )
      * \return  A newly created right vector
      */
    virtual Vector::shared_ptr getRightVector() const;

    /** \brief Get a left vector( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$\mathbf{y}\f$ is a left
     * vector )
      * \return  A newly created left vector
      */
    virtual Vector::shared_ptr getLeftVector() const;

    /** \brief Get the DOFManager associated with a right vector( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$\mathbf{x}\f$ is
     * a right vector )
      * \return  The DOFManager associated with a right vector
      */
    virtual Discretization::DOFManager::shared_ptr getRightDOFManager() const;

    /** \brief Get the DOFManager associated with a left vector( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$\mathbf{y}\f$ is
     * a left vector )
      * \return  The DOFManager associated with a left vector
      */
    virtual Discretization::DOFManager::shared_ptr getLeftDOFManager() const;

    /** \brief Compute the maximum column sum
      * \return  The L1 norm of the matrix
      */
    virtual double L1Norm() const;


protected:
    //! Unimplemented constructor
    DenseSerialMatrix();

    //! Protected copy constructor
    explicit DenseSerialMatrix( const DenseSerialMatrix & );

    /** \brief  Multiply two matrices and store in a third
     * \param[in]  other_op  The other matrix to multiply
     * \param[out] result  The matrix to store the result
     */
    virtual void multiply( shared_ptr other_op, shared_ptr &result );

    // AMP variables and DOFManagers for the left and right vectors
    AMP::LinearAlgebra::Variable::shared_ptr d_VariableLeft;
    AMP::LinearAlgebra::Variable::shared_ptr d_VariableRight;
    AMP::Discretization::DOFManager::shared_ptr d_DOFManagerLeft;
    AMP::Discretization::DOFManager::shared_ptr d_DOFManagerRight;

    // Data for the matrix
    size_t d_rows;
    size_t d_cols;
    double *d_M;
};
}
}


#endif
