#ifndef included_AMP_Vector
#define included_AMP_Vector


#include <iosfwd>
#include <string>

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/ParameterBase.h"
#include "AMP/utils/RNG.h"
#include "AMP/utils/enable_shared_from_this.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/data/VectorData.h"
#include "AMP/vectors/operations/VectorOperations.h"
#include <memory>


namespace AMP {
namespace LinearAlgebra {


class VectorSelector;


/** \brief Abstraction of a discrete Vector in a linear simulation
 * \details  This class encapsulates many BLAS level 1 operations
 * for use within AMP.  There are a variety of subclasses that
 * implement wrappers around various libraries so that these vectors
 * can be used by methods in those libraries.  Also, there are
 * several subclasses that allow combination of Vectors.  In general,
 * a user will only need to use Vector as the others are different
 * implementation details the developer need not manage.
 *
 * Vectors are created by factories.  For instance, SimpleVector
 * implements a create method that returns a shared pointer to a
 * Vector (Vector::shared_ptr).  The MeshAdpater interface implements
 * a createVector interface that also returns a Vector::shared_ptr.
 * Unless some specialized use of a Vector is needed (see SimpleVector),
 * the Vector::shared_ptr is the only type that should be used when
 * implementing math in AMP.
 *
 * Each Vector is associated with a Variable.  This Variable describes
 * the field this vector represents.
 * The problem may be discretized on a domain and linearized to
 * \f$\mathbf{L}\mathbf{\tilde{u}}=\mathbf{f}\f$.  In this case
 * \f$\mathbf{\tilde{u}}\f$ and \f$\mathbf{f}\f$ are Vectors.
 */
class Vector : virtual public VectorData, public AMP::enable_shared_from_this<Vector>
{

public: // typedefs
    using VectorData::ScatterType;
    using VectorData::UpdateState;

    /** \typedef shared_ptr
     * \brief Shorthand for shared pointer to Vector
     */
    typedef std::shared_ptr<Vector> shared_ptr;

    /** \typedef shared_ptr
     * \brief Shorthand for shared pointer to Vector
     */
    typedef std::shared_ptr<const Vector> const_shared_ptr;

    // Deprecated
    typedef VectorDataIterator<double> iterator;
    typedef VectorDataIterator<const double> const_iterator;


    // the next set of functions defines the public math. interface for vectors
public:
    /**
     * \brief  Set vector equal to x
     *      For Vectors, \f$\mathit{this}_i = x_i\f$.
     * \param[in] x         a vector
     */
    void copy( const VectorData &x );

    /**
     *\brief Set vector entries (including ghosts) to zero
     *\details This is equivalent (but more efficient) to calling setToScalar ( 0.0 ) followed by a
     *     makeConsistent(SET)
     */
    void zero( void );

    /**
     * \brief  Set all compenents of a vector to a scalar.
     *      For Vectors, the components of <em>this</em> are set to \f$\alpha\f$.
     * \param[in] alpha     scalar double value
     */
    void setToScalar( double alpha );

    /**
     * \brief Set data in this vector to random values on [0,1).
     */
    void setRandomValues( void );

    /**
     * \brief Set data in this vector to random values using
     *      a particular generator
     * \param[in] rng       The generator to use.
     */
    void setRandomValues( RNG::shared_ptr rng );

    /**
     * \brief  Set vector equal to scaled input.
     *      For Vectors, \f$\mathit{this}_i = \alpha x_i\f$.
     * \param[in] alpha     a scalar double
     * \param[in] x         a vector
     */
    void scale( double alpha, const VectorData &x );

    /**
     * \brief  Scale a vector.
     *     For Vectors, \f$\mathit{this}_i = \alpha\mathit{this}_i\f$.
     * \param[in] alpha     a scalar double
     */
    void scale( double alpha );

    /**
     * \brief  Adds two vectors.
     *      For Vectors, \f$\mathit{this}_i = x_i + y_i\f$.
     * \param[in] x         Input vector x
     * \param[in] y         Input vector y
     */
    void add( const VectorData &x, const VectorData &y );

    /**
     * \brief Subtracts one vector from another.
     *     For Vectors, \f$\mathit{this}_i = x_i - y_i\f$
     * \param[in] x         Input vector x
     * \param[in] y         Input vector y
     */
    void subtract( const VectorData &x, const VectorData &y );

    /**
     * \brief Component-wise multiply one vector with another.
     *    For Vectors, \f$\mathit{this}_i = x_i  y_i\f$
     * \param[in] x         Input vector x
     * \param[in] y         Input vector y
     */
    void multiply( const VectorData &x, const VectorData &y );

    /**
     * \brief Component-wise divide one vector by another.
     *    For Vectors, \f$\mathit{this}_i = x_i / y_i\f$
     * \param[in] x         Input vector x
     * \param[in] y         Input vector y
     */
    void divide( const VectorData &x, const VectorData &y );

    /**
     * \param x  a vector
     * \brief Set this to the component-wise reciprocal of a vector.  \f$\mathit{this}_i =
     * 1/x_i\f$.
     */
    void reciprocal( const VectorData &x );

    /**
     * \brief Set a vector to be a linear combination of two vectors.
     *      \f$\mathit{this}_i = \alpha x_i + \beta y_i\f$.
     * \param[in] alpha     a scalar
     * \param[in] x         a vector
     * \param[in] beta      a scalar
     * \param[in] y         a vector
     */
    void linearSum( double alpha, const VectorData &x, double beta, const VectorData &y );

    /**
     * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
     * \param[in] alpha    a scalar
     * \param[in] x        a vector
     * \param[in] y        a vector
     */
    void axpy( double alpha, const VectorData &x, const VectorData &y );

    /**
     * \brief Set this vector alpha * x + this.
     *     \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
     * \param[in] alpha    a scalar
     * \param[in] beta     a scalar
     * \param[in] x        a vector
     */
    void axpby( double alpha, double beta, const VectorData &x );

    /**
     * \brief Set this to the component-wise absolute value of a vector.
     *     \f$\mathit{this}_i = |x_i|\f$.
     * \param[in] x        a vector
     */
    void abs( const VectorData &x );

    /**
     * \brief set vector to \f$x + \alpha \bar{1}\f$.
     * \param[in] x a vector
     * \param[in] alpha a scalar
     * \details  for vectors, \f$\mathit{this}_i = x_i + \alpha\f$.
     */
    void addScalar( const VectorData &x, double alpha_in );

    /**
     * \brief Return the minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    double min( void ) const;

    /**
     * \brief Return the maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    double max( void ) const;

    /**
     * \brief Return discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    double L1Norm( void ) const;

    /**
     * \brief Return discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    double L2Norm( void ) const;

    /**
     * \brief Return the @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    double maxNorm( void ) const;
    /**
     * \brief Returns the minimum of the quotient of two vectors:
     *    \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     * \param[in] x a vector
     * \param[in] y a vector
     * \return \f[\min_{i,y_i\neq0} x_i/y_i\f]
     */
    double minQuotient( const VectorData &x ) const;

    /**
     * \brief Return a weighted norm of a vector
     * \param[in] x a vector
     * \param[in] y a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_i x^2_iy^2_i}{n}}\f]
     */
    double wrmsNorm( const VectorData &x, const VectorData &y ) const;

    /**
     * \brief Return a weighted norm of a subset of a vector
     * \param[in] x a vector
     * \param[in] y a vector
     * \param[in] mask a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_{i,\mathit{mask}_i>0} x^2_iy^2_i}{n}}\f]
     */
    double wrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y ) const;


    /**
     * \brief Return the dot product of this vector with the argument vector.
     * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     * \param[in] x        a vector
     */
    double dot( const VectorData &x ) const;

    bool equals( const VectorData &a, double tol ) const;


    /**
     * \brief Return the local minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    double localMin( void ) const;

    /**
     * \brief Return the local maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    double localMax( void ) const;

    /**
     * \brief Return local discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    double localL1Norm( void ) const;

    /**
     * \brief Return local discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    double localL2Norm( void ) const;

    /**
     * \brief Return the local @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    double localMaxNorm( void ) const;

    /**
     * \brief Return the local dot product of this vector with the argument vector.
     * \details Returns \f[\sum_i x_i \mathit{this}_i\f]
     * \param[in] x        a vector
     */
    double localDot( const VectorData &x ) const;

    /**
     * \brief Returns the local minimum of the quotient of two vectors:
     *    \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     * \param[in] x a vector
     * \return \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     */
    double localMinQuotient( const VectorData &x ) const;

    /**
     * \brief Return a weighted norm of a vector
     * \param[in] x a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_i x^2_i \mathit{this}^2_i}{n}}\f]
     */
    double localWrmsNorm( const VectorData &x ) const;

    /**
     * \brief Return a weighted norm of a subset of a vector
     * \param[in] x a vector
     * \param[in] mask a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_{i,\mathit{mask}_i>0}
     * \mathit{this}^2_iy^2_i}{n}}\f]
     */
    double
    localWrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y ) const;

    /**
     * \brief  Determine if the local portion of two vectors are equal using an absolute tolerance
     * \param[in] rhs      Vector to compare to
     * \param[in] tol      Tolerance of comparison
     * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
     */
    bool localEquals( const VectorData &x, double tol = 0.000001 ) const;


public: // shared_ptr wrappers
    /**
     * \brief  Determine if two vectors are equal using an absolute tolerance
     * \param[in] rhs      Vector to compare to
     * \param[in] tol      Tolerance of comparison
     * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
     */
    bool equals( std::shared_ptr<const VectorData> rhs, double tol = 0.000001 ) const;

    /// @copydoc VectorData::scale(double,const VectorData&)
    void scale( double alpha, std::shared_ptr<const VectorData> x );

    /// @copydoc VectorData::copy(const VectorData&)
    void copy( std::shared_ptr<const VectorData> x );

    /// @copydoc VectorData::add(const VectorData&,const VectorData&)
    void add( std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y );

    /// @copydoc VectorData::addScalar(const VectorData&,double)
    void addScalar( std::shared_ptr<const VectorData> x, double alpha );

    /// @copydoc VectorData::subtract(const VectorData&,const VectorData&)
    void subtract( std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y );

    /// @copydoc VectorData::multiply(const VectorData&,const VectorData&)
    void multiply( std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y );

    /// @copydoc VectorData::divide(const VectorData&,const VectorData&)
    void divide( std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y );

    /**
     * \param x  a vector
     * \brief Set this to the component-wise reciprocal of a vector.  \f$\mathit{this}_i =
     * 1/x_i\f$.
     */
    void reciprocal( std::shared_ptr<const VectorData> x );

    /**
     * \brief Set a vector to be a linear combination of two vectors.
     *      \f$\mathit{this}_i = \alpha x_i + \beta y_i\f$.
     * \param[in] alpha     a scalar
     * \param[in] x         a vector
     * \param[in] beta      a scalar
     * \param[in] y         a vector
     */
    void linearSum( double alpha,
                    std::shared_ptr<const VectorData> x,
                    double beta,
                    std::shared_ptr<const VectorData> y );

    /**
     * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
     * \param[in] alpha    a scalar
     * \param[in] x        a vector
     * \param[in] y        a vector
     */
    void
    axpy( double alpha, std::shared_ptr<const VectorData> x, std::shared_ptr<const VectorData> y );

    /**
     * \brief Set this vector alpha * x + this.
     *     \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
     * \param[in] alpha    a scalar
     * \param[in] beta     a scalar
     * \param[in] x        a vector
     */
    void axpby( double alpha, double beta, std::shared_ptr<const VectorData> x );

    /**
     * \brief Set this to the component-wise absolute value of a vector.
     *     \f$\mathit{this}_i = |x_i|\f$.
     * \param[in] x        a vector
     */
    void abs( std::shared_ptr<const VectorData> x );

    /**
     * \brief Return the dot product of this vector with the argument vector.
     * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     * \param[in] x        a vector
     */
    double dot( std::shared_ptr<const VectorData> x ) const;

    /**
     * \brief Returns the minimum of the quotient of two vectors:
     *    \f[\min_{i,y_i\neq0} x_i/\mathit{this}_i\f]
     * \param[in] x a vector
     * \param[in] y a vector
     * \return \f[\min_{i,y_i\neq0} x_i/y_i\f]
     */
    double minQuotient( std::shared_ptr<const VectorData> x ) const;

    /**
     * \brief Return a weighted norm of a vector
     * \param[in] x a vector
     * \param[in] y a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_i x^2_iy^2_i}{n}}\f]
     */
    double wrmsNorm( std::shared_ptr<const VectorData> x,
                     std::shared_ptr<const VectorData> y ) const;

    /**
     * \brief Return a weighted norm of a subset of a vector
     * \param[in] x a vector
     * \param[in] y a vector
     * \param[in] mask a vector
     * \return \f[\sqrt{\frac{\displaystyle \sum_{i,\mathit{mask}_i>0} x^2_iy^2_i}{n}}\f]
     */
    double wrmsNormMask( std::shared_ptr<const VectorData> x,
                         std::shared_ptr<const VectorData> mask,
                         std::shared_ptr<const VectorData> y ) const;

public: // Virtual functions
    /** \brief Return the name of the vector
     */
    virtual std::string type() const = 0;

    //! \name Vector memory manipulation
    //! \brief These methods control memory allocation, copying data, aliasing data, and swapping
    //! pointers among Vector
    //! instantiations
    //@{
    /** \brief Allocate space in the same fashion as <i>this</i>
     * \param[in] name  The variable to associate with the new vector
     * \details  This will allocate new space with identical layout as <i>this</i>.
     * \return  A Vector shared pointer
     * It will have the same number of blocks, each with the same engines and same number of
     * entries.
     */
    virtual Vector::shared_ptr cloneVector( const Variable::shared_ptr name ) const = 0;

    /** \brief  Swap the data in this Vector for another
      * \param[in]  other  Vector to swap data with
      * \details Effectively, this is
      * \code
      Vector *a;
      Vector *b;
      std::swap ( a, b );
        \endcode
      * without a and b exchanging pointers.
     */
    virtual void swapVectors( Vector &other ) = 0;

    /**  \brief  Make <i>this</i> be an alias of another vector
     *  \param[in]  other  Vector to be aliased
     *  \details  This will make <i>this</i> "point" to other.
     */
    virtual void aliasVector( Vector &other ) = 0;

    /**
     * \brief This method is used to implement the assemble interface
     * of PETSc.
     * \details  This method is empty except for instantiations of NativePetscVector
     */
    virtual void assemble() = 0;

    /** \brief Allocate space in the same fashion as <i>this</i>
     * \details  This will allocate new space with identical layout as <i>this</i>.
     * \return  A Vector shared pointer
     * It will have the same number of blocks, each with the same engines and same number of
     * entries.  The vector will
     * be associated with the same Variable.
     */
    virtual Vector::shared_ptr cloneVector() const;

    /** \brief Allocate space in the same fashion as <i>this</i>
     * \param[in] name  Name to give the variable associated with this vector
     * \details  This will allocate new space with identical layout as <i>this</i>.
     * \return  A Vector shared pointer
     * It will have the same number of blocks, each with the same engines and same
     * number of entries.  The vector will be associated with a clone of the same Variable with the
     * given name
     */
    virtual Vector::shared_ptr cloneVector( const std::string &name ) const;

    /** \brief Retrieve a sub-vector associated with a particular Variable
     * \param[in] name  Variable by which to retrieve a subvector
     * \return  A Vector shared pointer
     * \see MultiVector
     */
    virtual Vector::shared_ptr subsetVectorForVariable( Variable::const_shared_ptr name );

    /** \brief Retrieve a sub-vector associated with a particular Variable
     * \param[in] name  Variable by which to retrieve a subvector
     * \return  A Vector shared pointer
     * \see MultiVector
     */
    virtual Vector::const_shared_ptr
    constSubsetVectorForVariable( Variable::const_shared_ptr name ) const;

    /** \brief Return a parameters description of this vector
     * \return Parameters
     */
    virtual std::shared_ptr<ParameterBase> getParameters();

    /** \brief  Selects a portion of this vector and creates a view.
      * \param[in]  criterion  The method for deciding inclusion in the view
      * \param[in]  variable_name  The name of the vector to be created
      * \details To use, we recommend the following pattern
      \code
      // Vector to be "view"ed
      Vector::shared_ptr   data;

      // .. set up all the data storage in data

      // Get a view on the data tagged displacement
      Vector::shared_ptr  displacement = data->select ( VS_ByVariableName ( "displacement" ),
      "displacement view" );
      \endcode
      */
    virtual shared_ptr select( const VectorSelector &criterion, const std::string &variable_name );

    /** \brief  Selects a portion of this vector and creates a view.
      * \param[in]  criterion  The method for deciding inclusion in the view
      * \param[in]  variable_name  The name of the vector to be created
      * \details To use, we recommend the following pattern
      \code
      // Vector to be "view"ed
      Vector::shared_ptr   data;

      // .. set up all the data storage in data

      // Get a view on the data tagged displacement
      Vector::shared_ptr  displacement = data->select ( VS_ByVariableName ( "displacement" ),
      "displacement view" );
      \endcode
      */
    virtual const_shared_ptr constSelect( const VectorSelector &criterion,
                                          const std::string &variable_name ) const;

    //! Get the DOFManager for this Vector
    virtual AMP::Discretization::DOFManager::shared_ptr getDOFManager() const;

    /** \brief  Return the communicator this Vector spans
     */
    virtual AMP_MPI getComm() const override;


    /** \brief  Selects a portion of this vector and puts a view into a vector
     * \param[in]  criterion  The method for deciding inclusion in the view
     */
    virtual Vector::shared_ptr selectInto( const VectorSelector &criterion );

    // This is the const version of selectInto.
    virtual Vector::const_shared_ptr selectInto( const VectorSelector &criterion ) const;

    //! Return the pointer to the VectorData
    VectorData *getVectorData() { return d_VectorData; }

    //! Return the pointer to the VectorData
    const VectorData *getVectorData() const { return d_VectorData; }

    //! Return the pointer to the VectorOperation
    std::shared_ptr<VectorOperations> getVectorOperations() { return d_VectorOps; }

    //! Return the pointer to the VectorOperation
    std::shared_ptr<const VectorOperations> getVectorOperations() const { return d_VectorOps; }

    bool hasComm( void ) const override { return ( getVectorData()->getCommunicationList() != nullptr ); }

public: // Constructor/destructors
    /** \brief Constructor
     * \param[in] parameters  A pointer to a parameters class
     * \see VectorParameters
     */
    explicit Vector( VectorParameters::shared_ptr parameters );

    /** \brief Destructor
     */
    virtual ~Vector();


public: // Non-virtual functions
    void copyVector( std::shared_ptr<const Vector> x ) { d_VectorOps->copy( *x, *d_VectorData ); }

    /** \brief Change the variable associated with this vector
     * \param[in] name  The new variable
     */
    void setVariable( const Variable::shared_ptr name );

    /** \brief  Get the variable associated with this vector
     * \return  A shared point to the Variable associated with this Vector
     */
    const Variable::shared_ptr getVariable() const;

    /** \brief  Get the variable associated with this vector
     * \return  A shared point to the Variable associated with this Vector
     */
    Variable::shared_ptr getVariable();

    /** \brief Retrieve a sub-vector associated with a particular Variable
     * \param[in] name  Variable by which to retrieve a subvector
     * \return  A Vector shared pointer
     * \see MultiVector
     */
    Vector::shared_ptr subsetVectorForVariable( const std::string &name );


    /** \brief Retrieve a sub-vector associated with a particular Variable
     * \param[in] name  Variable by which to retrieve a subvector
     * \return  A Vector shared pointer
     * \see MultiVector
     */
    Vector::const_shared_ptr constSubsetVectorForVariable( const std::string &name ) const;

    /** \brief  Swap the data in this Vector for another
      * \param[in]  other Vector to swap data with
      * \details Effectively, this is
      * \code
      Vector *a;
      Vector *b;
      std::swap ( a, b );
        \endcode
      * without a and b exchanging pointers.
     */
    void swapVectors( Vector::shared_ptr other );

    /**  \brief  Make <i>this</i> be an alias of another vector
     *  \param[in]  other  Vector to be aliased
     *  \details  This will make <i>this</i> "point" to other.
     */
    void aliasVector( Vector::shared_ptr other );

    /** \brief  If a particular type of view of this Vector has been created,
     * return it.
     * \tparam VIEW_TYPE The type of view to look for
     * \return A view of this vector
     */
    template<typename VIEW_TYPE>
    Vector::shared_ptr getView() const;

    /** \brief  If a particular type of view of this Vector has been created,
     * return true.
     * \tparam VIEW_TYPE The type of view to look for
     * \return True if a view of this type has been created.  False, otherwise.
     */
    template<typename VIEW_TYPE>
    bool hasView() const;

    /** \brief Add a view of this vector to an internal queue.
     * \param[in] v The view to add
     */
    void registerView( Vector::shared_ptr v ) const;

    /** \brief Set the default RNG of this vector
     * \param[in] rng  The generator to set
     */
    static void setDefaultRNG( RNG::shared_ptr rng );

    /** \brief Get the current default RNG of this vector
     * \return  The current default RNG.
     * \details  If setDefaultRNG has not been called, this returns
     * an AMP::RNG base class.
     */
    static RNG::shared_ptr getDefaultRNG();

    /** \brief Associate the ghost buffer of a Vector with this Vector
     * \param in  The Vector to share a ghost buffer with
     */
    void aliasGhostBuffer( Vector::shared_ptr in );

protected: // Internal data
    // A default RNG to use when one is not specified
    static RNG::shared_ptr d_DefaultRNG;

    // The Variable associated with this Vector
    Variable::shared_ptr d_pVariable;

    // Constructor
    Vector();

    //! The DOF_Manager
    AMP::Discretization::DOFManager::shared_ptr d_DOFManager = nullptr;

    // Pointer to a VectorOperations object
    std::shared_ptr<VectorOperations> d_VectorOps = nullptr;

protected:
    // Pointer to *this as a VectorData object
    VectorData *d_VectorData = nullptr;
    
private:
    // The following are not implemented
    explicit Vector( const Vector & );
    void operator=( const Vector & );

    // output stream for vector data
    std::ostream *d_output_stream;

    std::shared_ptr<std::vector<std::weak_ptr<Vector>>> d_Views;
};


//! Stream operator
std::ostream &operator<<( std::ostream &out, const Vector::shared_ptr );
//! Stream operator
std::ostream &operator<<( std::ostream &out, const Vector & );


} // namespace LinearAlgebra
} // namespace AMP

#include "Vector.inline.h"

#endif
