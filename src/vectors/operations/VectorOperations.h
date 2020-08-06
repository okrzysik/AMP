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

    /**
     * \brief  Set vector equal to x
     *      For Vectors, \f$\mathit{this}_i = x_i\f$.
     * \param[in] x         a vector
     */
    virtual void copy( const VectorOperations &x ) = 0;

    /**
     *\brief Set vector entries (including ghosts) to zero
     *\details This is equivalent (but more efficient) to calling setToScalar ( 0.0 ) followed by a
     *     makeConsistent(SET)
     */
    virtual void zero() = 0;

    /**
     * \brief  Set all compenents of a vector to a scalar.
     *      For Vectors, the components of <em>this</em> are set to \f$\alpha\f$.
     * \param[in] alpha     scalar double value
     */
    virtual void setToScalar( double alpha ) = 0;

    /**
     * \brief Set data in this vector to random values on [0,1).
     */
    virtual void setRandomValues( void ) = 0;

    /**
     * \brief Set data in this vector to random values using
     *      a particular generator
     * \param[in] rng       The generator to use.
     */
    virtual void setRandomValues( RNG::shared_ptr rng ) = 0;

    /**
     * \brief  Set vector equal to scaled input.
     *      For Vectors, \f$\mathit{this}_i = \alpha x_i\f$.
     * \param[in] alpha     a scalar double
     * \param[in] x         a vector
     */
    virtual void scale( double alpha, const VectorOperations &x ) = 0;

    /**
     * \brief  Scale a vector.
     *     For Vectors, \f$\mathit{this}_i = \alpha\mathit{this}_i\f$.
     * \param[in] alpha     a scalar double
     */
    virtual void scale( double alpha ) = 0;

    /**
     * \brief set vector to \f$x + \alpha \bar{1}\f$.
     * \param[in] x a vector
     * \param[in] alpha a scalar
     * \details  for vectors, \f$\mathit{this}_i = x_i + \alpha\f$.
     */
    virtual void addScalar( const VectorOperations &x, double alpha ) = 0;

    /**
     * \brief  Adds two vectors.
     *      For Vectors, \f$\mathit{this}_i = x_i + y_i\f$.
     * \param[in] x         Input vector x
     * \param[in] y         Input vector y
     */
    virtual void add( const VectorOperations &x, const VectorOperations &y ) = 0;

    /**
     * \brief Subtracts one vector from another.
     *     For Vectors, \f$\mathit{this}_i = x_i - y_i\f$
     * \param[in] x         Input vector x
     * \param[in] y         Input vector y
     */
    virtual void subtract( const VectorOperations &x, const VectorOperations &y ) = 0;

    /**
     * \brief Component-wise multiply one vector with another.
     *    For Vectors, \f$\mathit{this}_i = x_i  y_i\f$
     * \param[in] x         Input vector x
     * \param[in] y         Input vector y
     */
    virtual void multiply( const VectorOperations &x, const VectorOperations &y ) = 0;

    /**
     * \brief Component-wise divide one vector by another.
     *    For Vectors, \f$\mathit{this}_i = x_i / y_i\f$
     * \param[in] x         Input vector x
     * \param[in] y         Input vector y
     */
    virtual void divide( const VectorOperations &x, const VectorOperations &y ) = 0;

    /**
     * \param x  a vector
     * \brief Set this to the component-wise reciprocal of a vector.  \f$\mathit{this}_i =
     * 1/x_i\f$.
     */
    virtual void reciprocal( const VectorOperations &x ) = 0;

    /**
     * \brief Set a vector to be a linear combination of two vectors.
     *      \f$\mathit{this}_i = \alpha x_i + \beta y_i\f$.
     * \param[in] alpha     a scalar
     * \param[in] x         a vector
     * \param[in] beta      a scalar
     * \param[in] y         a vector
     */
    virtual void linearSum( double alpha,
                            const VectorOperations &x,
                            double beta,
                            const VectorOperations &y ) = 0;

    /**
     * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
     * \param[in] alpha    a scalar
     * \param[in] x        a vector
     * \param[in] y        a vector
     */
    virtual void axpy( double alpha, const VectorOperations &x, const VectorOperations &y ) = 0;

    /**
     * \brief Set this vector alpha * x + this.
     *     \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
     * \param[in] alpha    a scalar
     * \param[in] beta     a scalar
     * \param[in] x        a vector
     */
    virtual void axpby( double alpha, double beta, const VectorOperations &x ) = 0;

    /**
     * \brief Set this to the component-wise absolute value of a vector.
     *     \f$\mathit{this}_i = |x_i|\f$.
     * \param[in] x        a vector
     */
    virtual void abs( const VectorOperations &x ) = 0;

    /**
     * \brief Return the minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    virtual double min( void ) const;

    /**
     * \brief Return the maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    virtual double max( void ) const;

    /**
     * \brief Return discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    virtual double L1Norm( void ) const;

    /**
     * \brief Return discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    virtual double L2Norm( void ) const;

    /**
     * \brief Return the @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    virtual double maxNorm( void ) const;

    /**
     * \brief Return the dot product of this vector with the argument vector.
     * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     * \param[in] x        a vector
     */
    virtual double dot( const VectorOperations &x ) const;

    /**
     * \brief Return the local minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    virtual double localMin( void ) const = 0;

    /**
     * \brief Return the local maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    virtual double localMax( void ) const = 0;

    /**
     * \brief Return local discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    virtual double localL1Norm( void ) const = 0;

    /**
     * \brief Return local discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    virtual double localL2Norm( void ) const = 0;

    /**
     * \brief Return the local @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    virtual double localMaxNorm( void ) const = 0;

    /**
     * \brief Return the local dot product of this vector with the argument vector.
     * \details Returns \f[\sum_i x_i \mathit{this}_i\f]
     * \param[in] x        a vector
     */
    virtual double localDot( const VectorOperations &x ) const = 0;

    /**
     * \brief  Determine if the local portion of two vectors are equal using an absolute tolerance
     * \param[in] rhs      Vector to compare to
     * \param[in] tol      Tolerance of comparison
     * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
     */
    virtual bool localEquals( const VectorOperations &rhs, double tol = 0.000001 ) const = 0;

    /**
     * \brief Returns the local minimum of the quotient of two vectors:
     *    \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     * \param[in] x a vector
     * \return \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     */
    virtual double localMinQuotient( const VectorOperations &x ) const = 0;

    /**
     * \brief Return a weighted norm of a vector
     * \param[in] x a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_i x^2_i \mathit{this}^2_i}{n}}\f]
     */
    virtual double localWrmsNorm( const VectorOperations &x ) const = 0;

    /**
     * \brief Return a weighted norm of a subset of a vector
     * \param[in] x a vector
     * \param[in] mask a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_{i,\mathit{mask}_i>0}
     * \mathit{this}^2_iy^2_i}{n}}\f]
     */
    virtual double localWrmsNormMask( const VectorOperations &x,
                                      const VectorOperations &mask ) const = 0;
public:
    // VectorData versions
    virtual void copy( const VectorData &x, VectorData &z ) = 0;
    virtual void zero( VectorData &x ) = 0;
    virtual void setToScalar( double alpha, VectorData &z ) = 0;
    virtual void setRandomValues( VectorData &x ) = 0;    
    virtual void setRandomValues( RNG::shared_ptr rng, VectorData &x ) = 0;    
    virtual void scale( double alpha, const VectorData &x, VectorData &y ) = 0;
    virtual void scale( double alpha, VectorData &x ) = 0;
    virtual void add( const VectorData &x, const VectorData &y, VectorData &z ) = 0;
    virtual void subtract( const VectorData &x, const VectorData &y, VectorData &z ) = 0;
    virtual void multiply( const VectorData &x, const VectorData &y, VectorData &z ) = 0;
    virtual void divide( const VectorData &x, const VectorData &y, VectorData &z ) = 0;
    virtual void reciprocal( const VectorData &x, VectorData &y ) = 0;
    virtual void linearSum( double alpha,
			   const VectorData &x,
			   double beta,
			   const VectorData &y,
			   VectorData &z) = 0;
    virtual void axpy( double alpha, const VectorData &x, const VectorData &y, VectorData &z ) = 0;
    virtual void axpby( double alpha, double beta, const VectorData &x, VectorData &y ) = 0;
    virtual void abs( const VectorData &x, VectorData &z ) = 0;
    virtual void addScalar( const VectorData &x, double alpha_in, VectorData &y ) = 0;

    virtual double min( const VectorData &x ) const;
    virtual double max( const VectorData &x ) const;
    virtual double L1Norm( const VectorData &x ) const;
    virtual double L2Norm( const VectorData &x  ) const;
    virtual double maxNorm( const VectorData &x ) const;
    virtual double dot( const VectorData &x, const VectorData &y ) const;
    virtual double localMin( const VectorData &x ) const = 0;
    virtual double localMax( const VectorData &x ) const = 0;
    virtual double localL1Norm( const VectorData &x ) const = 0;
    virtual double localL2Norm( const VectorData &x  ) const = 0;
    virtual double localMaxNorm( const VectorData &x ) const = 0;
    virtual double localDot( const VectorData &x, const VectorData &y ) const = 0;
    virtual double localMinQuotient( const VectorData &x, const VectorData &y ) const = 0;
    virtual double localWrmsNorm( const VectorData &x, const VectorData &y ) const = 0;
    virtual double localWrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y ) const = 0;
    virtual bool   localEquals( const VectorData &x, const VectorData &y, double tol = 0.000001 ) const = 0;
    
public: // Non-virtual functions
    /**
     * \brief  Determine if two vectors are equal using an absolute tolerance
     * \param[in] rhs      Vector to compare to
     * \param[in] tol      Tolerance of comparison
     * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
     */
    bool equals( const VectorOperations &rhs, double tol = 0.000001 ) const;


    /**
     * \brief Returns the minimum of the quotient of two vectors:
     *    \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     * \param[in] x a vector
     * \param[in] y a vector
     * \return \f[\min_{i,y_i\neq0} x_i/y_i\f]
     */
    static double minQuotient( const VectorOperations &x, const VectorOperations &y );

    /**
     * \brief Return a weighted norm of a vector
     * \param[in] x a vector
     * \param[in] y a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_i x^2_iy^2_i}{n}}\f]
     */
    static double wrmsNorm( const VectorOperations &x, const VectorOperations &y );

    /**
     * \brief Return a weighted norm of a subset of a vector
     * \param[in] x a vector
     * \param[in] y a vector
     * \param[in] mask a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_{i,\mathit{mask}_i>0} x^2_iy^2_i}{n}}\f]
     */
    static double wrmsNormMask( const VectorOperations &x,
                                const VectorOperations &y,
                                const VectorOperations &mask );


public: // shared_ptr wrappers
    /**
     * \brief  Determine if two vectors are equal using an absolute tolerance
     * \param[in] rhs      Vector to compare to
     * \param[in] tol      Tolerance of comparison
     * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
     */
    inline bool equals( std::shared_ptr<const VectorOperations> rhs, double tol = 0.000001 );
#if 1
    /// @copydoc VectorData::scale(double,const VectorData&)
    inline void scale( double alpha, std::shared_ptr<const VectorData> x, std::shared_ptr<VectorData> y );

    /// @copydoc VectorData::copy(const VectorData&)
    inline void copy( std::shared_ptr<const VectorData> x, std::shared_ptr<VectorData> y );

    /// @copydoc VectorData::add(const VectorData&,const VectorData&)
    inline void add( std::shared_ptr<const VectorData> x,
                     std::shared_ptr<const VectorData> y,
		     std::shared_ptr<VectorData> z );

    /// @copydoc VectorData::addScalar(const VectorData&,double)
    inline void addScalar( std::shared_ptr<const VectorData> x, double alpha, std::shared_ptr<VectorData> y );

    /// @copydoc VectorData::subtract(const VectorData&,const VectorData&)
    inline void subtract( std::shared_ptr<const VectorData> x,
                          std::shared_ptr<const VectorData> y,
			  std::shared_ptr<VectorData> z );

    /// @copydoc VectorData::multiply(const VectorData&,const VectorData&)
    inline void multiply( std::shared_ptr<const VectorData> x,
                          std::shared_ptr<const VectorData> y,
			  std::shared_ptr<VectorData> z );

    /// @copydoc VectorData::divide(const VectorData&,const VectorData&)
    inline void divide( std::shared_ptr<const VectorData> x,
                        std::shared_ptr<const VectorData> y,
			std::shared_ptr<VectorData> z );

    /**
     * \param x  a vector
     * \brief Set this to the component-wise reciprocal of a vector.  \f$\mathit{this}_i =
     * 1/x_i\f$.
     */
    inline void reciprocal( std::shared_ptr<const VectorData> x, std::shared_ptr<VectorData> y );

    /**
     * \brief Set a vector to be a linear combination of two vectors.
     *      \f$\mathit{this}_i = \alpha x_i + \beta y_i\f$.
     * \param[in] alpha     a scalar
     * \param[in] x         a vector
     * \param[in] beta      a scalar
     * \param[in] y         a vector
     */
    inline void linearSum( double alpha,
                           std::shared_ptr<const VectorData> x,
                           double beta,
                           std::shared_ptr<const VectorData> y,
			   std::shared_ptr<VectorData> z );

    /**
     * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
     * \param[in] alpha    a scalar
     * \param[in] x        a vector
     * \param[in] y        a vector
     */
    inline void axpy( double alpha,
                      std::shared_ptr<const VectorData> x,
                      std::shared_ptr<const VectorData> y,
		      std::shared_ptr<VectorData> z );

    /**
     * \brief Set this vector alpha * x + this.
     *     \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
     * \param[in] alpha    a scalar
     * \param[in] beta     a scalar
     * \param[in] x        a vector
     */
    inline void axpby( double alpha, double beta, std::shared_ptr<const VectorData> x, std::shared_ptr<VectorData> y );

    /**
     * \brief Set this to the component-wise absolute value of a vector.
     *     \f$\mathit{this}_i = |x_i|\f$.
     * \param[in] x        a vector
     */
    inline void abs( std::shared_ptr<const VectorData> x, std::shared_ptr<VectorData> y );

    /**
     * \brief Return the dot product of this vector with the argument vector.
     * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     * \param[in] x        a vector
     */
    inline double dot( std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y ) const;
#else

    /// @copydoc VectorOperations::scale(double,const VectorOperations&)
    inline void scale( double alpha, std::shared_ptr<const VectorOperations> x );

    /// @copydoc VectorOperations::copy(const VectorOperations&)
    inline void copy( std::shared_ptr<const VectorOperations> x );

    /// @copydoc VectorOperations::add(const VectorOperations&,const VectorOperations&)
    inline void add( std::shared_ptr<const VectorOperations> x,
                     std::shared_ptr<const VectorOperations> y );

    /// @copydoc VectorOperations::addScalar(const VectorOperations&,double)
    inline void addScalar( std::shared_ptr<const VectorOperations> x, double alpha );

    /// @copydoc VectorOperations::subtract(const VectorOperations&,const VectorOperations&)
    inline void subtract( std::shared_ptr<const VectorOperations> x,
                          std::shared_ptr<const VectorOperations> y );

    /// @copydoc VectorOperations::multiply(const VectorOperations&,const VectorOperations&)
    inline void multiply( std::shared_ptr<const VectorOperations> x,
                          std::shared_ptr<const VectorOperations> y );

    /// @copydoc VectorOperations::divide(const VectorOperations&,const VectorOperations&)
    inline void divide( std::shared_ptr<const VectorOperations> x,
                        std::shared_ptr<const VectorOperations> y );

    /**
     * \param x  a vector
     * \brief Set this to the component-wise reciprocal of a vector.  \f$\mathit{this}_i =
     * 1/x_i\f$.
     */
    inline void reciprocal( std::shared_ptr<const VectorOperations> x );

    /**
     * \brief Set a vector to be a linear combination of two vectors.
     *      \f$\mathit{this}_i = \alpha x_i + \beta y_i\f$.
     * \param[in] alpha     a scalar
     * \param[in] x         a vector
     * \param[in] beta      a scalar
     * \param[in] y         a vector
     */
    inline void linearSum( double alpha,
                           std::shared_ptr<const VectorOperations> x,
                           double beta,
                           std::shared_ptr<const VectorOperations> y );

    /**
     * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
     * \param[in] alpha    a scalar
     * \param[in] x        a vector
     * \param[in] y        a vector
     */
    inline void axpy( double alpha,
                      std::shared_ptr<const VectorOperations> x,
                      std::shared_ptr<const VectorOperations> y );

    /**
     * \brief Set this vector alpha * x + this.
     *     \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
     * \param[in] alpha    a scalar
     * \param[in] beta     a scalar
     * \param[in] x        a vector
     */
    inline void axpby( double alpha, double beta, std::shared_ptr<const VectorOperations> x );

    /**
     * \brief Set this to the component-wise absolute value of a vector.
     *     \f$\mathit{this}_i = |x_i|\f$.
     * \param[in] x        a vector
     */
    inline void abs( std::shared_ptr<const VectorOperations> x );

    /**
     * \brief Return the dot product of this vector with the argument vector.
     * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     * \param[in] x        a vector
     */
    inline double dot( std::shared_ptr<const VectorOperations> x ) const;

#endif
    /**
     * \brief Returns the minimum of the quotient of two vectors:
     *    \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     * \param[in] x a vector
     * \param[in] y a vector
     * \return \f[\min_{i,y_i\neq0} x_i/y_i\f]
     */
    static inline double minQuotient( std::shared_ptr<const VectorOperations> x,
                                      std::shared_ptr<const VectorOperations> y );

    /**
     * \brief Return a weighted norm of a vector
     * \param[in] x a vector
     * \param[in] y a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_i x^2_iy^2_i}{n}}\f]
     */
    static inline double wrmsNorm( std::shared_ptr<const VectorOperations> x,
                                   std::shared_ptr<const VectorOperations> y );

    /**
     * \brief Return a weighted norm of a subset of a vector
     * \param[in] x a vector
     * \param[in] y a vector
     * \param[in] mask a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_{i,\mathit{mask}_i>0} x^2_iy^2_i}{n}}\f]
     */
    static inline double wrmsNormMask( std::shared_ptr<const VectorOperations> x,
                                       std::shared_ptr<const VectorOperations> y,
                                       std::shared_ptr<const VectorOperations> mask );
public:
    //! Return the pointer to the VectorData
    inline VectorData *getVectorData() { return d_VectorData; }

    //! Return the pointer to the VectorData
    inline const VectorData *getVectorData() const { return d_VectorData; }

    //! Do we have a valid communicator
    inline bool hasComm() const;

    //! Do we have a valid communicator
    inline const AMP_MPI &getComm() const;

protected:
    VectorOperations();

    inline bool hasGhosts() const;
    inline std::vector<double> &getGhosts();


protected: // Internal data
    // Pointer to *this as a VectorData object
    VectorData *d_VectorData;
};


} // namespace LinearAlgebra
} // namespace AMP

#include "VectorOperations.inline.h"

#endif
