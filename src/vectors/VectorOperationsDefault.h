#ifndef included_AMP_VectorOperationsDefault
#define included_AMP_VectorOperationsDefault


#include "vectors/VectorOperations.h"


namespace AMP {
namespace LinearAlgebra {


/**
  * \brief  A default set of vector operations
  * \details VectorOperationsDefault impliments a default set of 
  *    vector operations on the CPU. 
  */
class VectorOperationsDefault : virtual public VectorOperations
{
public:

    // Constructor
    VectorOperationsDefault() {}

    //! Destructor
    virtual ~VectorOperationsDefault() {}


    /**
     * \param  alpha a scalar double
     * \brief  Set all compenents of a vector to a scalar.
     * For Vectors, the components of <em>this</em> are set to \f$\alpha\f$.
     */
    virtual void setToScalar( double alpha ) override = 0;

    /**
     * \param  alpha  a scalar double
     * \param  x  a vector
     * \brief  Set vector equal to scaled input.
     * For Vectors, \f$\mathit{this}_i = \alpha x_i\f$.
     */
    virtual void scale( double alpha, const VectorOperations &x ) override;

    /**
     * \param  alpha  a scalar double
     *
     * \brief  Scale a vector.
     * For Vectors, \f$\mathit{this}_i = \alpha\mathit{this}_i\f$.
     */
    virtual void scale( double alpha ) override;

    /**
     * \param  x  a vector
     * \param  y  a vector
     * \brief  Adds two vectors.
     * For Vectors, \f$\mathit{this}_i = x_i + y_i\f$.
     */
    virtual void add( const VectorOperations &x, const VectorOperations &y ) override;

    /**
      * \param x  a vector
      * \param y  a vector
      * \brief Subtracts one vector from another.
      * For Vectors, \f$\mathit{this}_i = x_i - y_i\f$
     */
    virtual void subtract( const VectorOperations &x, const VectorOperations &y ) override;

    /**
      * \param x  a vector
      * \param y  a vector
      * \brief Component-wise multiply one vector with another.
      * For Vectors, \f$\mathit{this}_i = x_i  y_i\f$
     */
    virtual void multiply( const VectorOperations &x, const VectorOperations &y ) override;

    /**
      * \param x  a vector
      * \param y  a vector
      * \brief Component-wise divide one vector by another.
      * For Vectors, \f$\mathit{this}_i = x_i / y_i\f$
     */
    virtual void divide( const VectorOperations &x, const VectorOperations &y ) override;

    /**
      * \param x  a vector
      * \brief Set this to the component-wise reciprocal of a vector.  \f$\mathit{this}_i =
     * 1/x_i\f$.
     */
    virtual void reciprocal( const VectorOperations &x ) override;


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
                            const VectorOperations &y ) override;

    /**
      * \param alpha a scalar
      * \param x a vector
      * \param y a vector
      * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
     */
    virtual void axpy( double alpha, const VectorOperations &x, const VectorOperations &y ) override;

    /**
      * \param alpha a scalar
      * \param beta a scalar
      * \param x  a vector
      * \brief Set this vector alpha * x + this.
      * \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
      */
    virtual void axpby( double alpha, double beta, const VectorOperations &x ) override;

    /**
      * \param x a vector
      * \brief Set this to the component-wise absolute value of a vector.
      * \f$\mathit{this}_i = |x_i|\f$.
     */
    virtual void abs( const VectorOperations &x ) override;

    /**
      * \brief Return the local minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    virtual double localMin( void ) const override;

    /**
      * \brief Return the local maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    virtual double localMax( void ) const override;

    /**
     * \brief Return local discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    virtual double localL1Norm( void ) const override;

    /**
     * \brief Return local discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    virtual double localL2Norm( void ) const override;

    /**
     * \brief Return the local @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    virtual double localMaxNorm( void ) const override;

    /**
      * \param[in] x a vector
      * \brief Return the local dot product of this vector with the argument vector.
      * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     */
    virtual double localDot( const VectorOperations& x ) const override;

    /**
     * \brief Set data in this vector to random values on [0,1).
     */
    virtual void setRandomValues( void ) override = 0;

};


} // LinearAlgebra namespace
} // AMP namespace

#endif
