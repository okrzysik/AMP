#ifndef included_AMP_VectorOperations
#define included_AMP_VectorOperations


#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/RNG.h"
#include <memory>
#include <vector>


namespace AMP {
namespace LinearAlgebra {


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
     * \param[in] alpha     scalar double value
     */
    virtual void setToScalar( double alpha, VectorData &z ) = 0;

    /**
     * \brief Set data in this vector to random values on [0,1).
     */
    virtual void setRandomValues( VectorData &x ) = 0;

    /**
     * \brief Set data in this vector to random values using
     *      a particular generator
     * \param[in] rng       The generator to use.
     */
    virtual void setRandomValues( RNG::shared_ptr rng, VectorData &x ) = 0;

    /**
     * \brief  Set vector equal to scaled input.
     *      For Vectors, \f$\mathit{this}_i = \alpha x_i\f$.
     * \param[in] alpha     a scalar double
     * \param[in] x         a vector
     */
    virtual void scale( double alpha, const VectorData &x, VectorData &y ) = 0;

    /**
     * \brief  Scale a vector.
     *     For Vectors, \f$\mathit{this}_i = \alpha\mathit{this}_i\f$.
     * \param[in] alpha     a scalar double
     */
    virtual void scale( double alpha, VectorData &x ) = 0;

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
    virtual void linearSum(
        double alpha, const VectorData &x, double beta, const VectorData &y, VectorData &z ) = 0;

    /**
     * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
     * \param[in] alpha    a scalar
     * \param[in] x        a vector
     * \param[in] y        a vector
     */
    virtual void axpy( double alpha, const VectorData &x, const VectorData &y, VectorData &z ) = 0;

    /**
     * \brief Set this vector alpha * x + this.
     *     \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
     * \param[in] alpha    a scalar
     * \param[in] beta     a scalar
     * \param[in] x        a vector
     */
    virtual void axpby( double alpha, double beta, const VectorData &x, VectorData &y ) = 0;

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
    virtual void addScalar( const VectorData &x, double alpha_in, VectorData &y ) = 0;

    /**
     * \brief Return the minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    virtual double min( const VectorData &x ) const;

    /**
     * \brief Return the maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    virtual double max( const VectorData &x ) const;

    /**
     * \brief Return discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    virtual double L1Norm( const VectorData &x ) const;

    /**
     * \brief Return discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    virtual double L2Norm( const VectorData &x ) const;

    /**
     * \brief Return the @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    virtual double maxNorm( const VectorData &x ) const;
    /**
     * \brief Returns the minimum of the quotient of two vectors:
     *    \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     * \param[in] x a vector
     * \param[in] y a vector
     * \return \f[\min_{i,y_i\neq0} x_i/y_i\f]
     */
    virtual double minQuotient( const VectorData &x, const VectorData &y ) const;

    /**
     * \brief Return a weighted norm of a vector
     * \param[in] x a vector
     * \param[in] y a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_i x^2_iy^2_i}{n}}\f]
     */
    virtual double wrmsNorm( const VectorData &x, const VectorData &y ) const;

    /**
     * \brief Return a weighted norm of a subset of a vector
     * \param[in] x a vector
     * \param[in] y a vector
     * \param[in] mask a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_{i,\mathit{mask}_i>0} x^2_iy^2_i}{n}}\f]
     */
    virtual double
    wrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y ) const;

    /**
     * \brief Return the dot product of this vector with the argument vector.
     * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     * \param[in] x        a vector
     */
    virtual double dot( const VectorData &x, const VectorData &y ) const;

    virtual bool equals( const VectorData &a, const VectorData &b, double tol ) const;

    /**
     * \brief Return the local minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    virtual double localMin( const VectorData &x ) const = 0;

    /**
     * \brief Return the local maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    virtual double localMax( const VectorData &x ) const = 0;

    /**
     * \brief Return local discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    virtual double localL1Norm( const VectorData &x ) const = 0;

    /**
     * \brief Return local discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    virtual double localL2Norm( const VectorData &x ) const = 0;

    /**
     * \brief Return the local @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    virtual double localMaxNorm( const VectorData &x ) const = 0;

    /**
     * \brief Return the local dot product of this vector with the argument vector.
     * \details Returns \f[\sum_i x_i \mathit{this}_i\f]
     * \param[in] x        a vector
     */
    virtual double localDot( const VectorData &x, const VectorData &y ) const = 0;

    /**
     * \brief Returns the local minimum of the quotient of two vectors:
     *    \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     * \param[in] x a vector
     * \return \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     */
    virtual double localMinQuotient( const VectorData &x, const VectorData &y ) const = 0;

    /**
     * \brief Return a weighted norm of a vector
     * \param[in] x a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_i x^2_i \mathit{this}^2_i}{n}}\f]
     */
    virtual double localWrmsNorm( const VectorData &x, const VectorData &y ) const = 0;

    /**
     * \brief Return a weighted norm of a subset of a vector
     * \param[in] x a vector
     * \param[in] mask a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_{i,\mathit{mask}_i>0}
     * \mathit{this}^2_iy^2_i}{n}}\f]
     */
    virtual double
    localWrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y ) const = 0;

    /**
     * \brief  Determine if the local portion of two vectors are equal using an absolute tolerance
     * \param[in] rhs      Vector to compare to
     * \param[in] tol      Tolerance of comparison
     * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
     */
    virtual bool
    localEquals( const VectorData &x, const VectorData &y, double tol = 0.000001 ) const = 0;


protected:
    VectorOperations();
};


} // namespace LinearAlgebra
} // namespace AMP

#endif
