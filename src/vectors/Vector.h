#ifndef included_AMP_Vector
#define included_AMP_Vector


#define DEPRECATED(x,y) std::cout << "Deprecated method: " << x << " use " << y << std::endl;

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#ifndef included_ParameterBase
#include "utils/ParameterBase.h"
#endif

#include "utils/PIO.h"
#include "boost/shared_ptr.hpp"
#include "boost/enable_shared_from_this.hpp"
#include <string>


#include "utils/RNG.h"
#include "utils/Castable.h"
#include "CommunicationList.h"

#include "Variable.h"
#include "VectorEngine.h"

#include "discretization/DOF_Manager.h"

#ifndef NULL
#define NULL (0)
#endif

namespace AMP {
namespace LinearAlgebra {


//! Parameters used to instantiate a Vector
class VectorParameters : public ParameterBase, public Castable
{
public:
    //! Convenience typedef
    typedef boost::shared_ptr<VectorParameters>   shared_ptr;

    //! The CommunicationList for a vector
    CommunicationList::shared_ptr             d_CommList;

    //! The DOF_Manager for a vector
    AMP::Discretization::DOFManager::shared_ptr   d_DOFManager;
};


class VectorDataIterator;
class ConstVectorDataIterator;
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

class Vector : virtual public VectorOperations, public boost::enable_shared_from_this<Vector>
{
public:
    /**\brief Flag to choose algorithm for makeConsistent
      *\see makeConsistent
      */
    enum ScatterType { CONSISTENT_ADD , CONSISTENT_SET };

    /** \typedef shared_ptr
      * \brief Shorthand for shared pointer to Vector
      */
    typedef boost::shared_ptr <Vector>     shared_ptr;

    /** \typedef shared_ptr
      * \brief Shorthand for shared pointer to Vector
      */
    typedef boost::shared_ptr <const Vector>     const_shared_ptr;

    /** \brief Return the name of the vector
      */
    virtual std::string type () const = 0;

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
    typedef VectorDataIterator         iterator;

    /** \typedef const_iterator
      * \brief An iterator for the data in a Vector
      */
    typedef ConstVectorDataIterator      const_iterator;
    /** \brief Constructor
      * \param[in] parameters  A pointer to a parameters class
      * \see VectorParameters
     */
    Vector( VectorParameters::shared_ptr  parameters);

    /** \brief Destructor
     */
    virtual ~Vector();

    /** \brief Change the variable associated with this vector
      * \param[in] name  The new variable
     */
    void setVariable(const Variable::shared_ptr name);

    /** \brief  Get the variable associated with this vector
      * \return  A shared point to the Variable associated with this Vector
      */
    const Variable::shared_ptr getVariable() const;

    /** \brief  Get the variable associated with this vector
      * \return  A shared point to the Variable associated with this Vector
      */
    Variable::shared_ptr getVariable();

    //! \name Vector memory manipulation
    //! \brief These methods control memory allocation, copying data, aliasing data, and swapping pointers among Vector instantiations
    //@{
    /** \brief Allocate space in the same fashion as <i>this</i>
      * \param[in] name  The variable to associate with the new vector
      * \details  This will allocate new space with identical layout as <i>this</i>.
      * \return  A Vector shared pointer
      * It will have the same number of blocks, each with the same engines and same number of entries.
     */
    virtual Vector::shared_ptr cloneVector(const Variable::shared_ptr name) const = 0;

    /** \brief Allocate space in the same fashion as <i>this</i>
      * \details  This will allocate new space with identical layout as <i>this</i>.
      * \return  A Vector shared pointer
      * It will have the same number of blocks, each with the same engines and same number of entries.  The vector will be associated with the same Variable.
     */
    virtual Vector::shared_ptr cloneVector () const;

    /** \brief Allocate space in the same fashion as <i>this</i>
      * \param[in] name  Name to give the variable associated with this vector
      * \details  This will allocate new space with identical layout as <i>this</i>.
      * \return  A Vector shared pointer
      * It will have the same number of blocks, each with the same engines and same 
      * number of entries.  The vector will be associated with a clone of the same Variable with the given name
     */
    virtual Vector::shared_ptr cloneVector ( const std::string &name ) const;

    /** \brief Retrieve a sub-vector associated with a particular Variable
      * \param[in] name  Variable by which to retrieve a subvector
      * \return  A Vector shared pointer
      * \see MultiVector
     */
    virtual Vector::shared_ptr subsetVectorForVariable ( const Variable::shared_ptr &name );

    /** \brief Obtain a particular contiguous block of data cast to RETURN_TYPE
      * \tparam RETURN_TYPE  The pointer type of the return
      * \param[in] i  Which block
      * \return A contiguous array of type RETURN_TYPE
      */
    template <typename RETURN_TYPE>
    RETURN_TYPE *  getRawDataBlock ( size_t i = 0 );

    /** \brief Obtain a particular contiguous block of data cast to RETURN_TYPE
      * \tparam RETURN_TYPE  The pointer type of the return
      * \param[in] i  Which block
      * \return A const contiguous array of type RETURN_TYPE
      */
    template <typename RETURN_TYPE>
    const RETURN_TYPE *  getRawDataBlock( size_t i = 0 ) const;

    /** \brief Number of blocks of contiguous data in the Vector
      * \return Number of blocks in the Vector
      * \details  A vector is not necessarily contiguous in memory.  This method
      * returns the number of contiguous blocks in memory used by this vector
      */
    virtual size_t  numberOfDataBlocks () const = 0;

    /** \brief Number of elements in a data block
      * \param[in] i  particular data block
      * \return The size of a particular block
      */
    virtual size_t  sizeOfDataBlock ( size_t i = 0 ) const = 0;

    /** \brief Copy the elements of a vector into <i>this</i>
      *   Note: if the ghosts in the rhs do not match the ghosts in this, 
      *   a makeConsistent is performed to fill the ghosts.  Otherwise it
      *   is assumed that rhs has consistent ghost values.
      * \param[in] rhs  a shared pointer to the Vector to copy the data from
     */
    virtual void copyVector ( const Vector::const_shared_ptr &rhs );

    /** \brief  Swap the data in this Bector for another
      * \param[in]  other  Vector to swap data with
      * \details Effectively, this is
      * \code
      Vector *a;
      Vector *b;
      std::swap ( a , b );
        \endcode
      * without a and b exchanging pointers.
     */
    virtual void swapVectors(Vector &other) = 0;

    /** \brief  Swap the data in this Vector for another
      * \param[in]  other Vector to swap data with
      * \details Effectively, this is
      * \code
      Vector *a;
      Vector *b;
      std::swap ( a , b );
        \endcode
      * without a and b exchanging pointers.
     */
    void swapVectors ( Vector::shared_ptr &other );

    /* \brief  Returns true if this vector has this element
     * \param[in]  GID  The global ID of the element
     */
    bool  containsGlobalElement ( size_t GID );


    /** \brief Return a parameters description of this vector
      * \return Parameters
      */
    boost::shared_ptr<ParameterBase> getParameters();

    /** \brief  Selects a portion of this vector and creates a view.
      * \param[in]  criterion  The method for deciding inclusion in the view
      * \param[in]  variable_name  The name of the vector to be created
      * \details To use, we recommend the following pattern
      \code
      // Vector to be "view"ed
      Vector::shared_ptr   data;

      // .. set up all the data storage in data

      // Get a view on the data tagged displacement
      Vector::shared_ptr  displacement = data->select ( VS_ByVariableName ( "displacement" ) , "displacement view" );
      \endcode
      */
    virtual shared_ptr   select ( const VectorSelector &criterion , const std::string &variable_name );

    /** \brief  Selects a portion of this vector and puts a view into a vector
      * \param[in]  criterion  The method for deciding inclusion in the view
      * \param[in,out]  vector  The vector to add the view to
      * \details  vector must be a MultiVector.  The easiest way to ensure this is to
      * create it with the select method.
      */
    virtual void       selectInto ( const VectorSelector &criterion , Vector::shared_ptr vector );


    /**  \brief  Make <i>this</i> be an alias of another vector
      *  \param[in]  other  Vector to be aliased
      *  \details  This will make <i>this</i> "point" to other.
     */
    virtual void aliasVector(Vector &other ) = 0;

    /**  \brief  Make <i>this</i> be an alias of another vector
      *  \param[in]  other  Vector to be aliased
      *  \details  This will make <i>this</i> "point" to other.
     */
    void aliasVector ( Vector::shared_ptr &other );

    /**\brief Copy data into this vector
      *\param[in] buf  Buffer to copy from
      */
    virtual void   putRawData ( double *buf ) = 0;

    /**\brief Copy data out of this vector
      *\param[out] buf  Buffer to copy to
      *\details The Vector will allocate *buf to the appropriate length
      */
    virtual void   copyOutRawData ( double **buf );

    /**\brief Number of elements "owned" by this core
      *\return  Number of entries stored contiguously on this processor
      *\details  For some types of variables, vectors may store "ghost"
      * data---possibly non-contiguous subsets of entries stored on other
      * cores.
      */
    virtual size_t getLocalSize() const = 0;

    /**\brief Number of total entries in this vector across all cores
      *\return Number of entries stored across all cores in this
      */
    virtual size_t getGlobalSize() const = 0;

    /**\brief The largest index in the vector (whether it is stored or not)
      *\details  Sparse vectors may not actually store the largest index
      * and getGlobalSize will only return the number of values stored
      *\return The largest index
      */
    virtual size_t getGlobalMaxID() const;

    /**\brief The largest index in the vector (whether it is stored or not)
      *\details  Sparse vectors may not actually store the largest index
      * and getGlobalSize will only return the number of values stored
      *\return The largest index
      */
    virtual size_t getLocalMaxID() const;

    virtual size_t getLocalStartID() const;

    /**\brief Number of entries "owned" by other cores stored on this
      * core.
      *\return Number of entries "owned" by other cores stored on this core
      */
    virtual size_t getGhostSize() const;
    //@}

    //! \name Linear Algebra Interface (Use this interface to perform math)
    //@{

    /**
     * \brief  Set all compenents of a vector to a scalar.
     * \param[in]  alpha a scalar double
     * For Vectors, the components of <em>this</em> are set to \f$\alpha\f$.
     */
    virtual void setToScalar(double alpha);

    /**
      *\brief Set vector entries (including ghosts) to zero
      *\details This is equivalent (but more efficient) to calling setToScalar ( 0.0 ) followed by a makeConsistent(SET)
      */
    void zero();

    /**
     * \brief  Scale a vector.
     * \param[in]  alpha  a scalar double
     *
     * For Vectors, \f$\mathit{this}_i = \alpha\mathit{this}_i\f$.
     */
    virtual void scale(double alpha);

    /**
     * \brief  Set vector equal to scaled input.
     * \param[in]  alpha  a scalar double
     * \param[in]  x  a vector
     * For Vectors, \f$\mathit{this}_i = \alpha x_i\f$.
     */
    void scale ( double alpha , const Vector::shared_ptr &x );

    /**
      * \brief set vector to \f$x + \alpha \bar{1}\f$.
      * \param[in] x a vector
      * \param[in] alpha a scalar
      * \details  for vectors, \f$\mathit{this}_i = x_i + \alpha\f$.
      */
    void addScalar ( const Vector::shared_ptr &x , double alpha );

    /**
     * \brief  Adds two vectors.
     * \param[in]  x  a vector
     * \param[in]  y  a vector
     * For Vectors, \f$\mathit{this}_i = x_i + y_i\f$.
     */
    void add ( const Vector::shared_ptr &x , const Vector::shared_ptr &y );

    /**
      * \brief Subtracts one vector from another.
      * \param[in] x  a vector
      * \param[in] y  a vector
      * For Vectors, \f$\mathit{this}_i = x_i - y_i\f$
     */
    void subtract ( const Vector::shared_ptr &x , const Vector::shared_ptr &y );

    /**
      * \brief Component-wise multiply one vector with another.
      * \param[in] x  a vector
      * \param[in] y  a vector
      * For Vectors, \f$\mathit{this}_i = x_i  y_i\f$
     */
    void multiply ( const Vector::shared_ptr &x , const Vector::shared_ptr &y );

    /**
      * \brief Component-wise divide one vector by another.
      * \param[in] x  a vector
      * \param[in] y  a vector
      * For Vectors, \f$\mathit{this}_i = x_i / y_i\f$
     */
    void divide ( const Vector::shared_ptr &x , const Vector::shared_ptr &y );

    /**
      * \brief Set this to the component-wise reciprocal of a vector.  \f$\mathit{this}_i = 1/x_i\f$.
      * \param[in] x  a vector
     */
    void reciprocal ( const Vector::shared_ptr &x );

    /**
     * \param[in] alpha a scalar
     * \param[in] x a vector
     * \param[in] beta a scalar
     * \param[in] y a vector
     * \brief Set a vector to be a linear combination of two vectors.
     *  \f$\mathit{this}_i = \alpha x_i + \beta y_i\f$.
     */
    void linearSum ( double alpha , const Vector::shared_ptr &x , double beta , const Vector::shared_ptr &y );

    /**
      * \param[in] alpha a scalar
      * \param[in] x a vector
      * \param[in] y a vector
      * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
     */
    void axpy ( double alpha , const Vector::shared_ptr &x , const Vector::shared_ptr &y );
    /**
      * \param[in] alpha a scalar
      * \param[in] beta a scalar
      * \param[in] x  a vector
      * \brief Set this vector alpha * x + this.
      * \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
      */
    void axpby( double alpha , double beta , const Vector::shared_ptr &x );

    /**
      * \param[in] x a vector
      * \brief Set this to the component-wise absolute value of a vector.
      * \f$\mathit{this}_i = |x_i|\f$.
     */
    void abs ( const Vector::shared_ptr &x );

    /**
     * \brief Set data in this vector to random values on [0,1).
     */
    virtual void setRandomValues(void);

    /**
     * \brief Set data in this vector to random values using
     * a particular generator
     * \param[in]  rng  The generator to use.
     */
    void  setRandomValues ( RNG::shared_ptr rng );

    /**
      * \brief Return the minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    virtual double min(void) const;

    /**
      * \brief Return the maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    virtual double max(void) const;

    /**
     * \brief Return discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    virtual double L1Norm(void) const;

    /**
     * \brief Return discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    virtual double L2Norm(void) const;

    /**
     * \brief Return the @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    virtual double maxNorm(void) const;

    /**
      * \param[in] x a vector
      * \brief Return the dot product of this vector with the argument vector.
      * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     */
    virtual double dot ( const Vector::shared_ptr &x );

    /**
      * \brief Return the local minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    virtual double localMin(void) const;

    /**
      * \brief Return the local maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    virtual double localMax(void) const;

    /**
     * \brief Return local discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    virtual double localL1Norm(void) const;

    /**
     * \brief Return local discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    virtual double localL2Norm(void) const;

    /**
     * \brief Return the local @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    virtual double localMaxNorm(void) const;

    /**
      * \param[in] x a vector
      * \brief Return the local dot product of this vector with the argument vector.
      * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     */
    virtual double localDot ( const boost::shared_ptr<const Vector> &x ) const;

    /**
      * \brief  Determine if two vectors are equal using an absolute tolerance
      * \param[in] rhs Vector to compare to
      * \param[in] tol Tolerance of comparison
      * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
      */
    bool  equals ( Vector::shared_ptr &rhs , double  tol = 0.000001 );

    /**
      * \fn equals (Vector & rhs , double tol )
      * \brief  Determine if two vectors are equal using an absolute tolerance
      * \param[in] rhs Vector to compare to
      * \param[in] tol Tolerance of comparison
      * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
      */
    virtual bool  equals ( Vector & rhs , double  tol = 0.000001 ); // Will be const one day soon
    //@}


    //! \name Manual vector manipulation
    //@{

    /**
      * \brief Return an iterator to the beginning of the data
      * \returns A VectorDataIterator
      * \details Since the Vector presents an interface to a contiguous
      * block of data, it is natural for it to provide a random access
      * iterator.
      * \warning If you are using PETSc and plan to write data to these
      * iterators, be sure to call DataChangeListener::dataChanged() on
      * the vector after use.
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
    virtual VectorDataIterator   begin ();
    /**
      * \brief Return an iterator to the beginning of the data
      * \returns A ConstVectorDataIterator
      * \details Since the Vector presents an interface to a contiguous
      * block of data, it is natural for it to provide a random access
      * iterator.
      * \warning If you are using PETSc and plan to write data to these
      * iterators, be sure to call DataChangeListener::dataChanged() on
      * the vector after use.
      */
    virtual ConstVectorDataIterator   begin () const;

    /**
      * \brief Return an iterator to the end of the data
      * \returns A VectorDataIterator
      * \details Since the Vector presents an interface to a contiguous
      * block of data, it is natural for it to provide a random access
      * iterator.
      * \warning If you are using PETSc and plan to write data to these
      * iterators, be sure to call DataChangeListener::dataChanged() on
      * the vector after use.
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
    virtual VectorDataIterator   end ();

    /**
      * \brief Return an iterator to the end of the data
      * \returns A ConstVectorDataIterator
      * \details Since the Vector presents an interface to a contiguous
      * block of data, it is natural for it to provide a random access
      * iterator.
      * \warning If you are using PETSc and plan to write data to these
      * iterators, be sure to call DataChangeListener::dataChanged() on
      * the vector after use.
      */
    virtual ConstVectorDataIterator   end () const;


    /**
      * \brief Set values in the vector by their local offset
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[in] vals the values to place in the vector
      * \details This will set the owned values for this core.  All indices are
      * from 0.
      * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
      */
    virtual void setValuesByLocalID ( int num , size_t *indices , const double *vals ) = 0;
    /**
      * \brief Set a single value in the vector by local ID
      * \param[in] i  offset of value to set
      * \param[in] val the value to place in the vector
      * \details An alias for setValuesByLocalID ( 1 , &num , &val );
      */
    void setValueByLocalID(size_t i, const double val);

    /**
      * \brief Set owned values using global identifier
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[in] vals the values to place in the vector
      *
      * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
      */
    virtual void setLocalValuesByGlobalID ( int num , size_t *indices , const double *vals ) = 0;

    /**
      * \brief Set a single owned value using global identifier
      * \param[in] i  offset of value to set
      * \param[in] val the value to place in the vector
      * \details An alias for setLocalValuesByGlobalID ( 1 , &i , &val );
      */
    void setLocalValueByGlobalID(size_t i, const double val);

    /**
      * \brief Set ghost values using global identifier
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[in] vals the values to place in the vector
      *
      * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
      */
    virtual void setGhostValuesByGlobalID ( int num , size_t *indices , const double *vals );

    /**
      * \brief Set a ghost owned value using global identifier
      * \param[in] i  offset of value to set
      * \param[in] val the value to place in the vector
      * \details An alias for setLocalValuesByGlobalID ( 1 , &i , &val );
      */
    void setGhostValueByGlobalID(size_t i, const double val);

    /**
      * \brief Set owned or shared values using global identifier
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[in] vals the values to place in the vector
      * \details Since the shared buffer and owned buffer are separate,
      * this function must sort the data by buffer before setting
      * values.
      */
    virtual void setValuesByGlobalID ( int num , size_t *indices , const double *vals );

    /**
      * \brief Set an owned or shared value using global identifier
      * \param[in] i  offset of value to set
      * \param[in] val the value to place in the vector
      * \details  An alias for setValuesByGlobalID ( 1 , &i , &val )
      */
    void setValueByGlobalID(size_t i, const double val);


    /**
      * \brief Add values to vector entities by their local offset
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[in] vals the values to place in the vector
      * \details This will set the owned values for this core.  All indices are
      * from 0.
      * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{this}_{\mathit{indices}_i} + \mathit{vals}_i \f$
      */
    virtual void addValuesByLocalID ( int num , size_t *indices , const double *vals ) = 0;
    /**
      * \brief Add a single value in the vector by local ID
      * \param[in] i  offset of value to set
      * \param[in] val the value to place in the vector
      * \details An alias for addValuesByLocalID ( 1 , &num , &val );
      */
    void addValueByLocalID(size_t i, const double val);

    /**
      * \brief Add owned values using global identifier
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[in] vals the values to place in the vector
      *
      * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{this}_{\mathit{indices}_i} + \mathit{vals}_i \f$
      */
    virtual void addLocalValuesByGlobalID ( int num , size_t *indices , const double *vals ) = 0;
    /**
      * \brief Add a single owned value using global identifier
      * \param[in] i  offset of value to set
      * \param[in] val the value to place in the vector
      * \details An alias for addLocalValuesByGlobalID ( 1 , &i , &val );
      */
    void addLocalValueByGlobalID(size_t i, const double val);

    /**
      * \brief Add owned or shared values using global identifier
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[in] vals the values to place in the vector
      * \details Since the shared buffer and owned buffer are separate,
      * this function must sort the data by buffer before setting
      * values.
      */
    virtual void addValuesByGlobalID ( int num , size_t *indices , const double *vals );

    /**
      * \brief Add an owned or shared value using global identifier
      * \param[in] i  offset of value to set
      * \param[in] val the value to place in the vector
      * \details  An alias for setValuesByGlobalID ( 1 , &i , &val )
      */
    void addValueByGlobalID(size_t i, const double val);

    /**
      * \brief get ghosted values to add to off-proc elements
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[in] vals the values to place in the vector
      * \details This will get the ghosted updates this processor has made.  All indices are
      * from global 0.
      */
    virtual void getGhostAddValuesByGlobalID ( int num , size_t *indices , double *vals ) const;

    /**
      * \brief get values in the vector by their local offset
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[out] vals the values to place in the vector
      * \details This will get the owned values for this core.  All indices are
      * from 0.
      */
    virtual void getValuesByGlobalID ( int num , size_t *indices , double *vals ) const;

    /**
      * \brief Return a value from the vector.
      * \param[in] i The global index into the vector
      * \return The value stored at the index
      * \details This uses getValuesByGlobalID to get the value
      */
    double getValueByGlobalID ( size_t i ) const;

    /**
      * \brief Get local values in the vector by their global offset
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[out] vals the values to place in the vector
      * \details This will get any value owned by this core.
      */
    virtual void getLocalValuesByGlobalID ( int num , size_t *indices , double *vals ) const = 0;

    /**
      * \brief Return a local value from the vector.
      * \param[in] i The global index into the vector
      * \return The value stored at the index
      * \details This uses getLocalValuesByGlobalID to get the value
      */
    double getLocalValueByGlobalID ( size_t i ) const;

    /**
      * \brief Get ghost values in the vector by their global offset
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[out] vals the values to place in the vector
      * \details This will get any value owned by this core.
      */
    virtual void getGhostValuesByGlobalID ( int num , size_t *indices , double *vals ) const;

    /**
      * \brief Return a ghost value from the vector.
      * \param[in] i The global index into the vector
      * \return The value stored at the index
      * \details This uses getGhostValuesByGlobalID to get the value
      */
    double getGhostValueByGlobalID ( size_t i ) const;


    /**
      * \brief Get local values in the vector by their global offset
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[out] vals the values to place in the vector
      * \details This will get any value used by this core.
      */
    virtual void getValuesByLocalID ( int num , size_t *indices , double *vals ) const;

    /**
      * \brief Return a local value from the vector.
      * \param[in] i The global index into the vector
      * \return The value stored at the index
      * \details This uses getValuesByGlobalID to get the value
      */
    double getValueByLocalID ( size_t i ) const;

    /**
      * \brief This method is used to implement the assemble interface
      * of PETSc.
      * \details  This method is empty except for instantiations of NativePetscVector
      */
    virtual void assemble() = 0;
    //@}



    //! \name VectorOperations virtual interface
    //@{
    virtual void scale(double alpha, const VectorOperations &x);
    virtual void add(const VectorOperations &x, const VectorOperations &y);
    virtual void subtract(const VectorOperations &x, const VectorOperations &y);
    virtual void multiply( const VectorOperations &x, const VectorOperations &y);
    virtual void divide( const VectorOperations &x, const VectorOperations &y);
    virtual void reciprocal(const VectorOperations &x);
    virtual void linearSum(double alpha, const VectorOperations &x,
          double beta, const VectorOperations &y);
    virtual void axpy(double alpha, const VectorOperations &x, const VectorOperations &y);
    virtual void axpby(double alpha, double beta, const VectorOperations &x);
    virtual void abs(const VectorOperations &x);
    virtual double dot(const VectorOperations &x) const;
    //@}

    //! \name Static methods for computation
    //@{
    /**
      * \brief Return a weighted norm of a vector
      * \param[in] x a vector
      * \param[in] y a vector
      * \return \f[\sqrt{\frac{\displaystyle \sum_i x^2_iy^2_i}{n}}\f]
      */
    static double wrmsNorm(const VectorOperations &x, const VectorOperations &y);

    /**
      * \brief Return a weighted norm of a subset of a vector
      * \param[in] x a vector
      * \param[in] y a vector
      * \param[in] mask a vector
      * \return \f[\sqrt{\frac{\displaystyle \sum_{i,\mathit{mask}_i>0} x^2_iy^2_i}{n}}\f]
      */

    static double wrmsNormMask(const VectorOperations &x, const VectorOperations &y, const VectorOperations &mask );

    /**
      * \brief Returns the minimum of the quotient of two vectors
      * \param[in] x a vector
      * \param[in] y a vector
      * \return \f[\min_{i,y_i\neq0} x_i/y_i\f]
      */

	    static double minQuotient(const VectorOperations &x, const VectorOperations &y);
    /**
      * \brief Returns the minimum of the quotient of two vectors
      * \param[in] x a vector
      * \param[in] y a vector
      * \return \f[\min_{i,y_i\neq0} x_i/y_i\f]
      */
    static double minQuotient(const Vector::shared_ptr &x, const Vector::shared_ptr &y);

    /**
      * \brief Return a weighted norm of a vector
      * \param[in] x a vector
      * \param[in] y a vector
      * \details Returns \f[\sqrt{\frac{\displaystyle \sum_i x^2_iy^2_i}{n}}\f]
      */
    static double wrmsNorm(const Vector::shared_ptr &x, const Vector::shared_ptr &y);


    //@}


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
    virtual void makeConsistent ( ScatterType  t );

    //! Get the CommunicationList for this Vector
    virtual CommunicationList::shared_ptr  getCommunicationList () const;

    //! Get the DOFManager for this Vector
    virtual AMP::Discretization::DOFManager::shared_ptr  getDOFManager () const;

    /**\brief Set the CommunicationList for this Vector
      *\details  Setting the CommunicationList for a Vector may involve
      * reallocating ghost storage.
      */
    virtual void setCommunicationList ( CommunicationList::shared_ptr  comm );

    /** \brief  Return the communicator this Vector spans
      */
    virtual AMP_MPI getComm() const;

    //@}


    /** \brief  If a particular type of view of this Vector has been created,
      * return it.
      * \tparam VIEW_TYPE The type of view to look for
      * \return A view of this vector
      */
    template <typename VIEW_TYPE>
    Vector::shared_ptr  getView () const;

    /** \brief  If a particular type of view of this Vector has been created,
      * return true.
      * \tparam VIEW_TYPE The type of view to look for
      * \return True if a view of this type has been created.  False, otherwise.
      */
    template <typename VIEW_TYPE>
    bool  hasView () const;

    /** \brief Add a view of this vector to an internal queue.
      * \param[in] v The view to add
      */
    void  registerView ( Vector::shared_ptr v ) const;

    /** \brief Write owned data to an std::ostream
      * \param[in] out  The output stream to write to.
      * \param[in] GIDoffset  A number to add to the global ID when writing information
      * \param[in] LIDoffset  A number to add to the local ID when writing information
      */
    virtual void  dumpOwnedData ( std::ostream &out , size_t GIDoffset=0 , size_t LIDoffset = 0 ) const;
    /** \brief Write data owned by other processors to an std::ostream
      * \param[in] out  The output stream to write to.
      * \param[in] offset  A number to add to the global ID when writing information
      */
    virtual void  dumpGhostedData ( std::ostream &out , size_t offset=0 ) const;

    /** \brief Set the default RNG of this vector
      * \param[in] rng  The generator to set
      */
    static  void  setDefaultRNG ( RNG::shared_ptr rng );

    /** \brief Get the current default RNG of this vector
      * \return  The current default RNG.
      * \details  If setDefaultRNG has not been called, this returns
      * an AMP::RNG base class.
      */
    RNG::shared_ptr  getDefaultRNG ();


    /**\brief The four states a Vector can be in
      *\see makeConsistent
      */
    enum UpdateState { UNCHANGED, LOCAL_CHANGED, ADDING, SETTING, MIXED };


    /** \brief  Return the current update state of the Vector
      * \details  This returns the effective update state of the 
      *  vector, including any vectors it contains.  The effective 
      *  state is defined as:
      *  UNCHANGED - All data and sub vectors are unchanged
      *  LOCAL_CHANGED - Local data may be modified, sub vectors must either
      *             be UNCHANGED or LOCAL_CHANGED.
      *  ADDING - Local and ghost data may be modified through add opperations, 
      *             sub vectors must be UNCHANGED, LOCAL_CHANGED, or ADDING
      *  SETTING - Local and ghost data may be modified through set opperations, 
      *             sub vectors must be UNCHANGED, LOCAL_CHANGED, or SETTING
      * If different subvectors have incompatible states ADDING and SETTING,
      * this function will return MIXED
      */
    virtual UpdateState  getUpdateStatus() const;


    /** \brief  Return the current update state of this Vector
      * \details  This returns the pointer to the update state
      *  of the current vector only (not vectors it contains).  
      *  It should NOT be used by users.
      */
    boost::shared_ptr<UpdateState>  getUpdateStatusPtr() const;


    /** \brief  Tie the current update state to another
      * \details  This sets the pointer to the update state
      *  of the current vector only (not vectors it contains).  
      *  It should NOT be used by users.
      * \param  rhs Pointer to share update state with
      */
    void setUpdateStatusPtr ( boost::shared_ptr<UpdateState> rhs );


protected:

    /** \brief  A default RNG to use when one is not specified
      */
    static   RNG::shared_ptr  d_DefaultRNG;


    /** \brief Copy ghosted vlues to a vector
      * \param[in] rhs  Vector to copy ghost values from
      * \details  This ensures that ghosted values on this and rhs
      * are the same without a call to makeConsistent.
      * \see makeConsistent
      */
    void  copyGhostValues ( const boost::shared_ptr<const Vector> &rhs );

    /** \brief Notify listeners that data has changed in this vector.
      */
    virtual void  dataChanged();

    /** \brief returns if two vectors are the same length globally and locally, otherwise throws an exception
      * \param rhs  Vector to compare with
      */
    void  requireSameSize ( Vector &rhs );

    /** \brief Return a pointer to a particular block of memory in the
      * vector
      * \param i The block to return
      */
    virtual void *getRawDataBlockAsVoid ( size_t i ) = 0;

    /** \brief Return a pointer to a particular block of memory in the
      * vector
      * \param i The block to return
      */
    virtual const void *getRawDataBlockAsVoid ( size_t i ) const = 0;

    /** \brief Helper function that probably doesn't need to exist anymore
      * \param  commList  The CommunicationList to add to the parameters
      * \details For internal use only
      */
    virtual void addCommunicationListToParameters ( CommunicationList::shared_ptr commList );

    /** \brief Associate the ghost buffer of a Vector with this Vector
      * \param in  The Vector to share a ghost buffer with
      */
    void  aliasGhostBuffer ( Vector::shared_ptr in );

    /**\brief  The communication list for this vector
      */
    CommunicationList::shared_ptr   d_CommList;

    /**\brief  The current update state for a vector
      *\details A Vector can be in one of three states. This is the current state of the vector
      * Because a vector can be composed of vectors, the update state needs to be shared between them
      */
    boost::shared_ptr<UpdateState>  d_UpdateState;

    /**\brief  The Variable associated with this Vector
      */
    Variable::shared_ptr        d_pVariable;

    /** \brief Constructor
    */
    Vector();

    //! The DOF_Manager
    AMP::Discretization::DOFManager::shared_ptr  d_DOFManager;

private:

    // The following are not implemented
    Vector(const Vector&);
    void operator=(const Vector&);

    // output stream for vector data
    std::ostream* d_output_stream;

    boost::shared_ptr<std::vector<double> >         d_Ghosts;
    boost::shared_ptr<std::vector<double> >         d_AddBuffer;
    boost::shared_ptr<std::vector<boost::weak_ptr<Vector> > >  d_Views;

};


//! Stream operator
std::ostream &operator << ( std::ostream &out , const Vector::shared_ptr );
//! Stream operator
std::ostream &operator << ( std::ostream &out , const Vector & );


}
}

#include "VectorDataIterator.h"
#include "Vector.inline.h"
#include "Vector.tmpl.h"
#endif


