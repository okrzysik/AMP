#ifndef included_AMP_VectorOperations
#define included_AMP_VectorOperations


#include "utils/shared_ptr.h"
#include <vector>

#include "utils/Castable.h"

namespace AMP {
namespace LinearAlgebra {

/**
  \brief  A temporary class used to disassociate data storage and mathematics
  from each other

  \details

  VectorOperations is a temporary class that is helping disassociate data storage
  and vector operations such as dot product, norms, etc.  Currently, there are
  two classes that inherit from VectorOperations:  Vector and VectorEngine.  At
  some time in the (not so) distant future, this class will be dissolved entirely
  as the VectorEngine class and the Vector class will have two distinct interfaces.
  Until then, the methods below will have two meanings, one for a Vector and one
  for a VectorEngine.

  Perhaps a word or two on the difference.  A Vector has data and a VectorEngine.
  A VectorEngine operates on data.  The difference can be seen in the Vec interface
  in PETSc.  A Vec holds data and keeps pointers to operation functions.  The
  engine is the litany of Vec functions:  VecAbs, VecSetValues, VecNorm, etc.

  If you are reading this portion of the documentation, odds are you do not need
  to know about VectorOperations.
  */

class VectorOperations : public Castable
{
public:
    /**
     * \brief  Virtual destructor
     *
     * Even though there are no members of this class, we declare a virtual destructor
     */
    virtual ~VectorOperations() {}


    /**
     * \param  alpha a scalar double
     * \brief  Set all compenents of a vector to a scalar.
     * For Vectors, the components of <em>this</em> are set to \f$\alpha\f$.
     */
    virtual void setToScalar( double alpha ) = 0;

    /**
     * \param  alpha  a scalar double
     * \param  x  a vector
     * \brief  Set vector equal to scaled input.
     * For Vectors, \f$\mathit{this}_i = \alpha x_i\f$.
     */
    virtual void scale( double alpha, const VectorOperations &x ) = 0;

    /**
     * \param  alpha  a scalar double
     *
     * \brief  Scale a vector.
     * For Vectors, \f$\mathit{this}_i = \alpha\mathit{this}_i\f$.
     */
    virtual void scale( double alpha ) = 0;

    /**
     * \param  x  a vector
     * \param  y  a vector
     * \brief  Adds two vectors.
     * For Vectors, \f$\mathit{this}_i = x_i + y_i\f$.
     */
    virtual void add( const VectorOperations &x, const VectorOperations &y ) = 0;

    /**
      * \param x  a vector
      * \param y  a vector
      * \brief Subtracts one vector from another.
      * For Vectors, \f$\mathit{this}_i = x_i - y_i\f$
     */
    virtual void subtract( const VectorOperations &x, const VectorOperations &y ) = 0;

    /**
      * \param x  a vector
      * \param y  a vector
      * \brief Component-wise multiply one vector with another.
      * For Vectors, \f$\mathit{this}_i = x_i  y_i\f$
     */
    virtual void multiply( const VectorOperations &x, const VectorOperations &y ) = 0;

    /**
      * \param x  a vector
      * \param y  a vector
      * \brief Component-wise divide one vector by another.
      * For Vectors, \f$\mathit{this}_i = x_i / y_i\f$
     */
    virtual void divide( const VectorOperations &x, const VectorOperations &y ) = 0;

    /**
      * \param x  a vector
      * \brief Set this to the component-wise reciprocal of a vector.  \f$\mathit{this}_i =
     * 1/x_i\f$.
     */
    virtual void reciprocal( const VectorOperations &x ) = 0;


    /**
     * \param alpha a scalar
     * \param x a vector
     * \param beta a scalar
     * \param y a vector
     * \brief Set a vector to be a linear combination of two vectors.
     *  \f$\mathit{this}_i = \alpha x_i + \beta y_i\f$.
     */
    virtual void linearSum( double alpha,
                            const VectorOperations &x,
                            double beta,
                            const VectorOperations &y ) = 0;

    /**
      * \param alpha a scalar
      * \param x a vector
      * \param y a vector
      * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
     */
    virtual void axpy( double alpha, const VectorOperations &x, const VectorOperations &y ) = 0;

    /**
      * \param alpha a scalar
      * \param beta a scalar
      * \param x  a vector
      * \brief Set this vector alpha * x + this.
      * \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
      */
    virtual void axpby( double alpha, double beta, const VectorOperations &x ) = 0;

    /**
      * \param x a vector
      * \brief Set this to the component-wise absolute value of a vector.
      * \f$\mathit{this}_i = |x_i|\f$.
     */
    virtual void abs( const VectorOperations &x ) = 0;

    /**
      * \brief Return the minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    virtual double min( void ) const = 0;

    /**
      * \brief Return the maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    virtual double max( void ) const = 0;

    /**
     * \brief Set data in this vector to random values on [0,1).
     */
    virtual void setRandomValues( void ) = 0;

    /**
     * \brief Return discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    virtual double L1Norm( void ) const = 0;

    /**
     * \brief Return discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    virtual double L2Norm( void ) const = 0;

    /**
     * \brief Return the @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    virtual double maxNorm( void ) const = 0;

    /**
      * \param x a vector
      * \brief Return the dot product of this vector with the argument vector.
      * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     */
    virtual double dot( const VectorOperations &x ) const = 0;
};
}
}

#endif
