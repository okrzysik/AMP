#ifndef included_AMP_Matrix
#define included_AMP_Matrix

#include "utils/Castable.h"
#include "utils/ParameterBase.h"
#include "vectors/Vector.h"
#include "Matrix.h"  //Does Matrix.h need to include Matrix.h?

#include "boost/shared_ptr.hpp"

namespace AMP {
namespace LinearAlgebra {

  /** \class MatrixParameters
    * \brief  Description of a matrix for construciton
    */
  class MatrixParameters : public Castable
  {
    public:
      /** \brief Convenience typedef
        */
      typedef boost::shared_ptr<MatrixParameters>     shared_ptr;
  };

  /** \class Matrix
    * \brief  An abstract interface for using and manipulating
    * matrices
    * \details  There are several different varieties of distributed
    * memory matrices.  While most operations between the varieties
    * can be abstracted away from the user, some cannot.  For this
    * reason, most of the time, this class will suffice as the
    * way to interact with a matrix.  Matrix creation may require
    * use of one of the derived classes.
    */
  class Matrix : public Castable 
  {
    public:
      //! Convenience typedef
      typedef  boost::shared_ptr<Matrix>      shared_ptr;

    protected:
      /** \brief Unimplemented constructor
        */
      Matrix();


      /** \brief Unused copy constructor
        */
      Matrix ( const Matrix & );

      /** \brief  Multiply two matrices and store in a third
        * \param[in]  other_op  The other matrix to multiply
        * \param[out] result  The matrix to store the result
        */
      virtual void multiply ( shared_ptr other_op , shared_ptr &result ) = 0;

    public:

      /** \brief Constructor
        * \param[in] params  Description of the matrix
        */
      Matrix( MatrixParameters::shared_ptr  params );

      /** \brief Destructor
        */
      virtual ~Matrix();

      /** \brief  Matrix-vector multiplication
        * \param[in]  in  The vector to multiply
        * \param[out] out The resulting vectory
        * \details  Compute \f$\mathbf{Ain} = \mathbf{out}\f$.
        */
      virtual void mult(const boost::shared_ptr<Vector> & in, 
          boost::shared_ptr<Vector> & out) = 0;

      /** \brief  Matrix transpose-vector multiplication
        * \param[in]  in  The vector to multiply
        * \param[out] out The resulting vectory
        * \details  Compute \f$\mathbf{A}^T\mathbf{in} = \mathbf{out}\f$.
        */
      virtual void multTranspose(const boost::shared_ptr<Vector> & in, 
          boost::shared_ptr<Vector> & out) = 0;


      /** \brief  Return a new matrix that is the transpose of this one
        * \return  A copy of this matrix transposed.
        */
      virtual shared_ptr  transpose () const;

      /** \brief  Return a matrix with the same non-zero and distributed structure
        * \return  The new matrix
        */
      virtual shared_ptr  cloneMatrix () const = 0;

      /** \brief  Scale the matrix by a scalar
        * \param[in] alpha  The value to scale by
        * \details  Compute \f$\mathbf{A} = \alpha\mathbf{A}\f$
        */
      virtual void  scale ( double alpha ) = 0;


      /** \brief  Compute the product of two matrices
        * \param[in] A  A multiplicand
        * \param[in] B  A multiplicand
        * \return The product \f$\mathbf{AB}\f$.
        */
      static shared_ptr  matMultiply ( shared_ptr A , shared_ptr B );

      /** \brief  Compute the linear combination of two matrices
        * \param[in] alpha  scalar
        * \param[in] X matrix
        * \details  Compute \f$\mathbf{THIS} = \alpha\mathbf{X} + \mathbf{THIS}\f$
        */
      virtual void  axpy ( double alpha , const Matrix &X ) = 0;

      /** \brief  Compute the linear combination of two matrices
        * \param[in] alpha  scalar
        * \param[in] X matrix
        * \details  Compute \f$\mathbf{THIS} = \alpha\mathbf{X} + \mathbf{THIS}\f$
        */
      void  axpy ( double alpha , const Matrix::shared_ptr  &X );


      /** \brief  Add values to those in the matrix
        * \param[in] num_rows The number of rows represented in values
        * \param[in] num_cols The number of cols represented in values
        * \param[in] rows  The row ids of values
        * \param[in] cols  The column ids of values
        * \param[in] values  The values to add to the matrix
        * \details  This method may fail if the matrix has not
        * allocated a particular (row,col) specified, depending
        * on the actual subclass of matrix used.
        */
      virtual void  addValuesByGlobalID ( int   num_rows ,
                                          int   num_cols ,
                                          int  *rows ,
                                          int  *cols ,
                                          double  *values ) = 0;

      /** \brief  Set values in the matrix
        * \param[in] num_rows The number of rows represented in values
        * \param[in] num_cols The number of cols represented in values
        * \param[in] rows  The row ids of values
        * \param[in] cols  The column ids of values
        * \param[in] values  The values to set to the matrix
        * \details  This method may fail if the matrix has not
        * allocated a particular (row,col) specified, depending
        * on the actual subclass of matrix used.
        */
      virtual void  setValuesByGlobalID ( int   num_rows ,
                                          int   num_cols ,
                                          int  *rows ,
                                          int  *cols ,
                                          double  *values ) = 0;


      /** \brief  Add values to those in the matrix
        * \param[in] row  The row id of value
        * \param[in] col  The column id of value
        * \param[in] value  The value to add to the matrix
        * \details  This method may fail if the matrix has not
        * allocated a particular (row,col) specified, depending
        * on the actual subclass of matrix used.
        */
      virtual void addValueByGlobalID ( int row , int col , double value );

      /** \brief  Set values in the matrix
        * \param[in] row  The row id of value
        * \param[in] col  The column id of value
        * \param[in] value  The value to set to the matrix
        * \details  This method may fail if the matrix has not
        * allocated a particular (row,col) specified, depending
        * on the actual subclass of matrix used.
        */
      virtual void setValueByGlobalID ( int row , int col , double value );

      /** \brief  Set the non-zeros of the matrix to a scalar
        * \param[in]  alpha  The value to set the non-zeros to
        */
      virtual void setScalar ( double alpha ) = 0;

      /** \brief  Set the non-zeros of the matrix to zero
        * \details  May not deallocate space.
        */
      void         zero ();

      /** \brief  Retrieve a row of the matrix in compressed format
        * \param[in]  row Which row
        * \param[out] cols  The column ids of the returned values
        * \param[out] values  The values in the row
        */
      virtual void getRowByGlobalID ( int row , std::vector<unsigned int> &cols , 
                                                  std::vector<double> &values ) const = 0;

      /** \brief  Set the diagonal to the values in a vector
        * \param[in] in The values to set the diagonal to
        */
      virtual void setDiagonal ( const Vector::shared_ptr &in ) = 0;

      /** \brief  Perform communication to ensure values in the
        * matrix are the same across cores.
        */
      virtual void makeConsistent () = 0;

      /** \brief  Get the number of rows in the matrix
        * \return  The number of rows
        */
      virtual size_t numRows () = 0;

      /** \brief  Get the number of columns in the matrix
        * \return  The number of columns
        */
      virtual size_t numColumns () = 0;

      /** \brief  Extract the diagonal from a matrix
        * \param[in]  buf  An optional vector to use as a buffer
        * \return  A vector of the diagonal values
        */
      virtual Vector::shared_ptr  extractDiagonal ( Vector::shared_ptr buf = Vector::shared_ptr() ) = 0;

      /** \brief Get a right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$\mathbf{x}\f$ is a right vector )
        * \return  A newly created right vector
        */
      virtual Vector::shared_ptr  getRightVector () = 0;

      /** \brief Get a left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$\mathbf{y}\f$ is a left vector )
        * \return  A newly created left vector
        */
      virtual Vector::shared_ptr  getLeftVector () = 0;

      /** \brief Compute the maximum column sum
        * \return  The L1 norm of the matrix
        */
      virtual double  L1Norm() const = 0;
  };

}
}

#include "Matrix.inline.h"
#endif


