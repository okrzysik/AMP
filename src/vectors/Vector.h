#ifndef included_AMP_Vector
#define included_AMP_Vector


#include <iostream>
#include <string>

#include "discretization/DOF_Manager.h"
#include "utils/Castable.h"
#include "utils/ParameterBase.h"
#include "utils/RNG.h"
#include "utils/enable_shared_from_this.h"
#include "utils/shared_ptr.h"
#include "vectors/Variable.h"
#include "vectors/VectorData.h"
#include "vectors/VectorOperations.h"
#include "vectors/VectorOperationsDefault.h"

namespace AMP {
namespace LinearAlgebra {


//! Parameters used to instantiate a Vector
class VectorParameters : public ParameterBase, public Castable
{
public:
    //! Convenience typedef
    typedef AMP::shared_ptr<VectorParameters> shared_ptr;

    //! The CommunicationList for a vector
    CommunicationList::shared_ptr d_CommList;

    //! The DOF_Manager for a vector
    AMP::Discretization::DOFManager::shared_ptr d_DOFManager;
};


class VectorDataIterator;
class ConstVectorDataIterator;
class VectorSelector;
class MultiVector;
class ManagedVector;


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
class Vector :
    virtual public VectorData,
    virtual public VectorOperations,
    virtual public VectorOperationsDefault,
    public AMP::enable_shared_from_this<Vector>
{

public: // typedefs

    using VectorData::ScatterType;
    using VectorData::UpdateState;

    /** \typedef shared_ptr
      * \brief Shorthand for shared pointer to Vector
      */
    typedef AMP::shared_ptr<Vector> shared_ptr;

    /** \typedef shared_ptr
      * \brief Shorthand for shared pointer to Vector
      */
    typedef AMP::shared_ptr<const Vector> const_shared_ptr;

    /** \typedef iterator
      * \brief An iterator for the data in a Vector---DO NOT USE WITHOUT READING DETAILS.
      * \see DataChangeFirer
      * \warning  If you understand what DataChangeFirer does and why, then you can use
      * the non-const iterator.
      * This can be used with the following pattern:
      * \code
      void square ( Vector::shared_ptr  vector )
      {
        Vector::iterator    cur_entry = vector->begin();
        while ( cur_entry != vector->end() )
        {
        (*cur_entry) = (*cur_entry)*(*cur_entry);
        cur_entry++;
        }
        if ( vector->isA<DataChangeFirer>() )
        {
        vector->castTo<DataChangeFirer>().fireDataChange();
        }
      }
        \endcode
      */
    typedef VectorDataIterator iterator;

    /** \typedef const_iterator
      * \brief An iterator for the data in a Vector
      */
    typedef ConstVectorDataIterator const_iterator;


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
    virtual Vector::shared_ptr subsetVectorForVariable( const Variable::shared_ptr &name );

    /** \brief Retrieve a sub-vector associated with a particular Variable
      * \param[in] name  Variable by which to retrieve a subvector
      * \return  A Vector shared pointer
      * \see MultiVector
     */
    virtual Vector::const_shared_ptr
    constSubsetVectorForVariable( const Variable::shared_ptr &name ) const;

    /** \brief Copy the elements of a vector into <i>this</i>
      *   Note: if the ghosts in the rhs do not match the ghosts in this,
      *   a makeConsistent is performed to fill the ghosts.  Otherwise it
      *   is assumed that rhs has consistent ghost values.
      * \param[in] rhs  a shared pointer to the Vector to copy the data from
     */
    virtual void copyVector( Vector::const_shared_ptr rhs );

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


    /** \brief Return a parameters description of this vector
      * \return Parameters
      */
    virtual AMP::shared_ptr<ParameterBase> getParameters();

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


    //! \name Parallel Management
    //@{

    /**
      * \brief Update shared values on entire communicator
      * \param t The type of scatter used to compute values
      * \details  There are two algorithms used by makeConsistent
      * - If t = CONSISTENT_SET, then owned values are
      *   sent to processors that share the value.  Shared values are
      *   overwritten
      * - If t = CONSISTENT_ADD, then shared values are accumulated
      *   on the core that owns it and applied as determined, either
      *   add or set.  Then, the values are broadcast out.
      *
      * Generally, when adding to a vector, the GATHER_SCATTER should
      * be used to make consistent.  When setting entries in a vector
      * the BROADCAST should be used.
      */
    virtual void makeConsistent( ScatterType t );


    //! Get the DOFManager for this Vector
    virtual AMP::Discretization::DOFManager::shared_ptr getDOFManager() const;

    /**\brief Set the CommunicationList for this Vector
      *\details  Setting the CommunicationList for a Vector may involve
      * reallocating ghost storage.
      */
    virtual void setCommunicationList( CommunicationList::shared_ptr comm );

    /** \brief  Return the communicator this Vector spans
      */
    virtual AMP_MPI getComm() const;



    /**
      * \fn equals (Vector & const rhs, double tol )
      * \brief  Determine if two vectors are equal using an absolute tolerance
      * \param[in] rhs Vector to compare to
      * \param[in] tol Tolerance of comparison
      * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
      */
    virtual bool equals( Vector const &rhs, double tol = 0.000001 ) const;


    /**
     * \brief Set data in this vector to random values on [0,1).
     */
    virtual void setRandomValues( void ) override;


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

    /**
     * \brief Set data in this vector to random values using
     * a particular generator
     * \param[in]  rng  The generator to use.
     */
    void setRandomValues( RNG::shared_ptr rng );


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
    inline Vector::shared_ptr subsetVectorForVariable( const std::string& name );


    /** \brief Retrieve a sub-vector associated with a particular Variable
      * \param[in] name  Variable by which to retrieve a subvector
      * \return  A Vector shared pointer
      * \see MultiVector
     */
    inline Vector::const_shared_ptr
    constSubsetVectorForVariable( const std::string& name ) const;

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

    /**
      *\brief Set vector entries (including ghosts) to zero
      *\details This is equivalent (but more efficient) to calling setToScalar ( 0.0 ) followed by a
      *makeConsistent(SET)
      */
    void zero();

    /**
     * \brief  Set vector equal to scaled input.
     * \param[in]  alpha  a scalar double
     * \param[in]  x  a vector
     * For Vectors, \f$\mathit{this}_i = \alpha x_i\f$.
     */
    void scale( double alpha, Vector::const_shared_ptr x );

    /**
      * \brief set vector to \f$x + \alpha \bar{1}\f$.
      * \param[in] x a vector
      * \param[in] alpha a scalar
      * \details  for vectors, \f$\mathit{this}_i = x_i + \alpha\f$.
      */
    void addScalar( Vector::const_shared_ptr x, double alpha );

    /**
     * \brief  Adds two vectors.
     * \param[in]  x  a vector
     * \param[in]  y  a vector
     * For Vectors, \f$\mathit{this}_i = x_i + y_i\f$.
     */
    void add( Vector::const_shared_ptr x, Vector::const_shared_ptr y );

    /**
      * \brief Subtracts one vector from another.
      * \param[in] x  a vector
      * \param[in] y  a vector
      * For Vectors, \f$\mathit{this}_i = x_i - y_i\f$
     */
    void subtract( Vector::const_shared_ptr x, Vector::const_shared_ptr y );

    /**
      * \brief Component-wise multiply one vector with another.
      * \param[in] x  a vector
      * \param[in] y  a vector
      * For Vectors, \f$\mathit{this}_i = x_i  y_i\f$
     */
    void multiply( Vector::const_shared_ptr x, Vector::const_shared_ptr y );

    /**
      * \brief Component-wise divide one vector by another.
      * \param[in] x  a vector
      * \param[in] y  a vector
      * For Vectors, \f$\mathit{this}_i = x_i / y_i\f$
     */
    void divide( Vector::const_shared_ptr x, Vector::const_shared_ptr y );

    /**
      * \brief Set this to the component-wise reciprocal of a vector.  \f$\mathit{this}_i =
     * 1/x_i\f$.
      * \param[in] x  a vector
     */
    void reciprocal( Vector::const_shared_ptr x );

    /**
     * \param[in] alpha a scalar
     * \param[in] x a vector
     * \param[in] beta a scalar
     * \param[in] y a vector
     * \brief Set a vector to be a linear combination of two vectors.
     *  \f$\mathit{this}_i = \alpha x_i + \beta y_i\f$.
     */
    void
    linearSum( double alpha, Vector::const_shared_ptr x, double beta, Vector::const_shared_ptr y );


    /**
      * \param[in] alpha a scalar
      * \param[in] x a vector
      * \param[in] y a vector
      * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
      *    Note: after this call, the data may not be in a consistent state,
      *    and may require calling makeConsistent(SET) if consistency is required.
     */
    void axpy( double alpha, Vector::const_shared_ptr x, Vector::const_shared_ptr y );

    /**
      * \param[in] alpha a scalar
      * \param[in] beta a scalar
      * \param[in] x  a vector
      * \brief Set this vector alpha * x + this.
      *    \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
      *    Note: after this call, the data may not be in a consistent state,
      *    and may require calling makeConsistent(SET) if consistency is required.
      */
    void axpby( double alpha, double beta, Vector::const_shared_ptr x );

    /**
      * \param[in] x a vector
      * \brief Set this to the component-wise absolute value of a vector.
      * \f$\mathit{this}_i = |x_i|\f$.
     */
    void abs( Vector::const_shared_ptr x );

    /**
      * \brief  Determine if two vectors are equal using an absolute tolerance
      * \param[in] rhs Vector to compare to
      * \param[in] tol Tolerance of comparison
      * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
      */
    bool equals( Vector::const_shared_ptr rhs, double tol = 0.000001 ) const;

    /**
      * \param[in] x a vector
      * \brief Return the dot product of this vector with the argument vector.
      * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     */
    inline double dot( Vector::const_shared_ptr x ) const;


    //! \name Static methods for computation
    //@{
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

    /**
      * \brief Returns the minimum of the quotient of two vectors
      * \param[in] x a vector
      * \param[in] y a vector
      * \return \f[\min_{i,y_i\neq0} x_i/y_i\f]
      */

    static double minQuotient( const VectorOperations &x, const VectorOperations &y );
    /**
      * \brief Returns the minimum of the quotient of two vectors
      * \param[in] x a vector
      * \param[in] y a vector
      * \return \f[\min_{i,y_i\neq0} x_i/y_i\f]
      */
    static double minQuotient( Vector::const_shared_ptr x, Vector::const_shared_ptr y );

    /**
      * \brief Return a weighted norm of a vector
      * \param[in] x a vector
      * \param[in] y a vector
      * \details Returns \f[\sqrt{\frac{\displaystyle \sum_i x^2_iy^2_i}{n}}\f]
      */
    static double wrmsNorm( Vector::const_shared_ptr x, Vector::const_shared_ptr y );

    /** \brief  If a particular type of view of this Vector has been created,
      * return it.
      * \tparam VIEW_TYPE The type of view to look for
      * \return A view of this vector
      */
    template <typename VIEW_TYPE>
    Vector::shared_ptr getView() const;

    /** \brief  If a particular type of view of this Vector has been created,
      * return true.
      * \tparam VIEW_TYPE The type of view to look for
      * \return True if a view of this type has been created.  False, otherwise.
      */
    template <typename VIEW_TYPE>
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
    RNG::shared_ptr getDefaultRNG();



protected: // Virtual functions

    /** \brief  Selects a portion of this vector and puts a view into a vector
      * \param[in]  criterion  The method for deciding inclusion in the view
      * \param[in,out]  vector  The vector to add the view to
      * \details  vector must be a MultiVector.  The easiest way to ensure this is to
      * create it with the select method.
      */
    virtual Vector::shared_ptr selectInto( const VectorSelector &criterion );

    // This is the const version of selectInto.
    virtual Vector::const_shared_ptr selectInto( const VectorSelector &criterion ) const;

    /** \brief Helper function that probably doesn't need to exist anymore
      * \param  commList  The CommunicationList to add to the parameters
      * \details For internal use only
      */
    virtual void addCommunicationListToParameters( CommunicationList::shared_ptr commList );


protected: // Non-virtual functions


    /** \brief Copy ghosted vlues to a vector
      * \param[in] rhs  Vector to copy ghost values from
      * \details  This ensures that ghosted values on this and rhs
      * are the same without a call to makeConsistent.
      * \see makeConsistent
      */
    void copyGhostValues( const AMP::shared_ptr<const Vector> &rhs );

    /** \brief returns if two vectors are the same length globally and locally, otherwise throws an
     * exception
      * \param rhs  Vector to compare with
      */
    void requireSameSize( Vector &rhs );

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
    AMP::Discretization::DOFManager::shared_ptr d_DOFManager;

    friend class ManagedVector;
    friend class MultiVector;

private:
    // The following are not implemented
    explicit Vector( const Vector & );
    void operator=( const Vector & );

    // output stream for vector data
    std::ostream *d_output_stream;

    AMP::shared_ptr<std::vector<AMP::weak_ptr<Vector>>> d_Views;


public: // Default implimentations (MOVE)

    virtual void setToScalar(double) override;
    virtual void dumpOwnedData( std::ostream &out, size_t GIDoffset=0, size_t LIDoffset=0 ) const override;
    virtual void dumpGhostedData( std::ostream &out, size_t offset=0 ) const override;


public: // Pull VectorOperations into the current scope
    using VectorOperationsDefault::add;
    using VectorOperationsDefault::abs;
    using VectorOperationsDefault::axpy;
    using VectorOperationsDefault::axpby;
    using VectorOperationsDefault::divide;
    using VectorOperationsDefault::dot;
    using VectorOperationsDefault::linearSum;
    using VectorOperationsDefault::multiply;
    using VectorOperationsDefault::scale;
    using VectorOperationsDefault::subtract;
    using VectorOperationsDefault::reciprocal;
};


//! Stream operator
std::ostream &operator<<( std::ostream &out, const Vector::shared_ptr );
//! Stream operator
std::ostream &operator<<( std::ostream &out, const Vector & );


} // LinearAlgebra namespace
} // AMP namespace

#include "Vector.inline.h"

#endif
