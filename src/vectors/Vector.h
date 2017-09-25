#ifndef included_AMP_Vector
#define included_AMP_Vector


#include <iostream>
#include <string>

#include "discretization/DOF_Manager.h"
#include "utils/ParameterBase.h"
#include "utils/RNG.h"
#include "utils/enable_shared_from_this.h"
#include "utils/shared_ptr.h"
#include "vectors/Variable.h"
#include "vectors/data/VectorData.h"
#include "vectors/operations/VectorOperations.h"


namespace AMP {
namespace LinearAlgebra {


//! Parameters used to instantiate a Vector
class VectorParameters : public ParameterBase
{
public:
    //! Convenience typedef
    typedef AMP::shared_ptr<VectorParameters> shared_ptr;

    //! The CommunicationList for a vector
    CommunicationList::shared_ptr d_CommList;

    //! The DOF_Manager for a vector
    AMP::Discretization::DOFManager::shared_ptr d_DOFManager;
};


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
class Vector :
    virtual public VectorData,
    virtual public VectorOperations,
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

    // Deprecated
    typedef VectorDataIterator<double> iterator;
    typedef VectorDataIterator<const double> const_iterator;


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
    virtual Vector::shared_ptr cloneVector( const std::string& name ) const;

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


    /** \brief  Selects a portion of this vector and puts a view into a vector
      * \param[in]  criterion  The method for deciding inclusion in the view
      */
    virtual Vector::shared_ptr selectInto( const VectorSelector &criterion );

    // This is the const version of selectInto.
    virtual Vector::const_shared_ptr selectInto( const VectorSelector &criterion ) const;

    /** \brief Helper function that probably doesn't need to exist anymore
      * \param  commList  The CommunicationList to add to the parameters
      * \details For internal use only
      */
    virtual void addCommunicationListToParameters( CommunicationList::shared_ptr commList );


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

    /// @copydoc VectorOperations::copy(const VectorOperations&)
    inline void copyVector( AMP::shared_ptr<const Vector> x ) { copy(*x); }

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
    AMP::Discretization::DOFManager::shared_ptr d_DOFManager;


private:
    // The following are not implemented
    explicit Vector( const Vector & );
    void operator=( const Vector & );

    // output stream for vector data
    std::ostream *d_output_stream;

    AMP::shared_ptr<std::vector<AMP::weak_ptr<Vector>>> d_Views;

};


//! Stream operator
std::ostream &operator<<( std::ostream &out, const Vector::shared_ptr );
//! Stream operator
std::ostream &operator<<( std::ostream &out, const Vector & );


} // LinearAlgebra namespace
} // AMP namespace

#include "Vector.inline.h"

#endif
