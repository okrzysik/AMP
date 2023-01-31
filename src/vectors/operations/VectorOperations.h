#ifndef included_AMP_VectorOperations
#define included_AMP_VectorOperations

#include "AMP/vectors/Scalar.h"

#include <memory>
#include <vector>


namespace AMP::LinearAlgebra {


class VectorData;


/**
 * \brief  A class used to hold vector operations
 * \details  VectorOperations is a class used to perform vector operations
 *    such as dot product, norms, etc.
 */
class VectorOperations
{
public:
    //! Destructor
    virtual ~VectorOperations() {}

    //! Get the type name
    virtual std::string VectorOpName() const = 0;

    //! Clone the operations
    virtual std::shared_ptr<VectorOperations> cloneOperations() const = 0;

    // VectorData versions
    /**
     * \brief  Set vector equal to x
     *      For Vectors, \f$\mathit{this}_i = x_i\f$.
     * \param[in] x         a vector
     */
    virtual void copy( const VectorData &x, VectorData &z ) = 0;

    /**
     *\brief Set vector entries (including ghosts) to zero
     *\details This is equivalent (but more efficient) to calling setToScalar ( 0.0 ) followed by a
     *     makeConsistent(SET)
     */
    virtual void zero( VectorData &x ) = 0;

    /**
     * \brief  Set all compenents of a vector to a scalar.
     *      For Vectors, the components of <em>this</em> are set to \f$\alpha\f$.
     * \param[in] alpha     scalar value
     */
    virtual void setToScalar( const Scalar &alpha, VectorData &z ) = 0;

    /**
     * \brief Set data in this vector to random values
     */
    virtual void setRandomValues( VectorData &x ) = 0;

    /**
     * \brief  Set vector equal to scaled input.
     *      For Vectors, \f$\mathit{this}_i = \alpha x_i\f$.
     * \param[in] alpha     a scalar
     * \param[in] x         a vector
     */
    virtual void scale( const Scalar &alpha, const VectorData &x, VectorData &y ) = 0;

    /**
     * \brief  Scale a vector.
     *     For Vectors, \f$\mathit{this}_i = \alpha\mathit{this}_i\f$.
     * \param[in] alpha     a scalar
     */
    virtual void scale( const Scalar &alpha, VectorData &x ) = 0;

    /**
     * \brief  Adds two vectors.
     *      For Vectors, \f$\mathit{this}_i = x_i + y_i\f$.
     * \param[in] x         Input vector x
     * \param[in] y         Input vector y
     */
    virtual void add( const VectorData &x, const VectorData &y, VectorData &z ) = 0;

    /**
     * \brief Subtracts one vector from another.
     *     For Vectors, \f$\mathit{this}_i = x_i - y_i\f$
     * \param[in] x         Input vector x
     * \param[in] y         Input vector y
     */
    virtual void subtract( const VectorData &x, const VectorData &y, VectorData &z ) = 0;

    /**
     * \brief Component-wise multiply one vector with another.
     *    For Vectors, \f$\mathit{this}_i = x_i  y_i\f$
     * \param[in] x         Input vector x
     * \param[in] y         Input vector y
     */
    virtual void multiply( const VectorData &x, const VectorData &y, VectorData &z ) = 0;

    /**
     * \brief Component-wise divide one vector by another.
     *    For Vectors, \f$\mathit{this}_i = x_i / y_i\f$
     * \param[in] x         Input vector x
     * \param[in] y         Input vector y
     */
    virtual void divide( const VectorData &x, const VectorData &y, VectorData &z ) = 0;

    /**
     * \param x  a vector
     * \brief Set this to the component-wise reciprocal of a vector.  \f$\mathit{this}_i =
     * 1/x_i\f$.
     */
    virtual void reciprocal( const VectorData &x, VectorData &y ) = 0;

    /**
     * \brief Set a vector to be a linear combination of two vectors.
     *      \f$\mathit{this}_i = \alpha x_i + \beta y_i\f$.
     * \param[in] alpha     a scalar
     * \param[in] x         a vector
     * \param[in] beta      a scalar
     * \param[in] y         a vector
     */
    virtual void linearSum( const Scalar &alpha,
                            const VectorData &x,
                            const Scalar &beta,
                            const VectorData &y,
                            VectorData &z ) = 0;

    /**
     * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
     * \param[in] alpha    a scalar
     * \param[in] x        a vector
     * \param[in] y        a vector
     */
    virtual void
    axpy( const Scalar &alpha, const VectorData &x, const VectorData &y, VectorData &z ) = 0;

    /**
     * \brief Set this vector alpha * x + this.
     *     \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
     * \param[in] alpha    a scalar
     * \param[in] beta     a scalar
     * \param[in] x        a vector
     */
    virtual void
    axpby( const Scalar &alpha, const Scalar &beta, const VectorData &x, VectorData &y ) = 0;

    /**
     * \brief Set this to the component-wise absolute value of a vector.
     *     \f$\mathit{this}_i = |x_i|\f$.
     * \param[in] x        a vector
     */
    virtual void abs( const VectorData &x, VectorData &z ) = 0;

    /**
     * \brief set vector to \f$x + \alpha \bar{1}\f$.
     * \param[in] x a vector
     * \param[in] alpha a scalar
     * \details  for vectors, \f$\mathit{this}_i = x_i + \alpha\f$.
     */
    virtual void addScalar( const VectorData &x, const Scalar &alpha_in, VectorData &y ) = 0;

    /**
     * \brief Return the minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    virtual Scalar min( const VectorData &x ) const;

    /**
     * \brief Return the maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    virtual Scalar max( const VectorData &x ) const;

    /**
     * \brief Return the sum of the values of the vector.
     */
    virtual Scalar sum( const VectorData &x ) const;

    /**
     * \brief Return the mean of the values of the vector.
     */
    Scalar mean( const VectorData &x ) const;

    /**
     * \brief Return discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    virtual Scalar L1Norm( const VectorData &x ) const;

    /**
     * \brief Return discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    virtual Scalar L2Norm( const VectorData &x ) const;

    /**
     * \brief Return the @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    virtual Scalar maxNorm( const VectorData &x ) const;
    /**
     * \brief Returns the minimum of the quotient of two vectors:
     *    \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     * \param[in] x a vector
     * \param[in] y a vector
     * \return \f[\min_{i,y_i\neq0} x_i/y_i\f]
     */
    virtual Scalar minQuotient( const VectorData &x, const VectorData &y ) const;

    /**
     * \brief Return a weighted norm of a vector
     * \param[in] x a vector
     * \param[in] y a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_i x^2_iy^2_i}{n}}\f]
     */
    virtual Scalar wrmsNorm( const VectorData &x, const VectorData &y ) const;

    /**
     * \brief Return a weighted norm of a subset of a vector
     * \param[in] x a vector
     * \param[in] y a vector
     * \param[in] mask a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_{i,\mathit{mask}_i>0} x^2_iy^2_i}{n}}\f]
     */
    virtual Scalar
    wrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y ) const;

    /**
     * \brief Return the dot product of this vector with the argument vector.
     * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     * \param[in] x        a vector
     */
    virtual Scalar dot( const VectorData &x, const VectorData &y ) const;

    /**
     * \brief Check if two vectors are equal
     * \details Returns true if all values are equal within tolerance
     * \param[in] x        a vector
     */
    virtual bool equals( const VectorData &a, const VectorData &b, const Scalar &tol ) const;

    /**
     * \brief Return the local minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    virtual Scalar localMin( const VectorData &x ) const = 0;

    /**
     * \brief Return the local maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    virtual Scalar localMax( const VectorData &x ) const = 0;

    /**
     * \brief Return the local sumof the vector.
     */
    virtual Scalar localSum( const VectorData &x ) const = 0;

    /**
     * \brief Return local discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    virtual Scalar localL1Norm( const VectorData &x ) const = 0;

    /**
     * \brief Return local discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    virtual Scalar localL2Norm( const VectorData &x ) const = 0;

    /**
     * \brief Return the local @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    virtual Scalar localMaxNorm( const VectorData &x ) const = 0;

    /**
     * \brief Return the local dot product of this vector with the argument vector.
     * \details Returns \f[\sum_i x_i \mathit{this}_i\f]
     * \param[in] x        a vector
     */
    virtual Scalar localDot( const VectorData &x, const VectorData &y ) const = 0;

    /**
     * \brief Returns the local minimum of the quotient of two vectors:
     *    \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     * \param[in] x a vector
     * \return \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     */
    virtual Scalar localMinQuotient( const VectorData &x, const VectorData &y ) const = 0;

    /**
     * \brief Return a weighted norm of a vector
     * \param[in] x a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_i x^2_i \mathit{this}^2_i}{n}}\f]
     */
    virtual Scalar localWrmsNorm( const VectorData &x, const VectorData &y ) const = 0;

    /**
     * \brief Return a weighted norm of a subset of a vector
     * \param[in] x a vector
     * \param[in] mask a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_{i,\mathit{mask}_i>0}
     * \mathit{this}^2_iy^2_i}{n}}\f]
     */
    virtual Scalar
    localWrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y ) const = 0;

    /**
     * \brief  Determine if the local portion of two vectors are equal using an absolute tolerance
     * \param[in] rhs      Vector to compare to
     * \param[in] tol      Tolerance of comparison
     * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
     */
    virtual bool
    localEquals( const VectorData &x, const VectorData &y, const Scalar &tol = 1e-6 ) const = 0;


    //! Get a unique id hash for the vector
    uint64_t getID() const;


protected:
    VectorOperations();
};


} // namespace AMP::LinearAlgebra

#endif
