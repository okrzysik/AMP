#ifndef included_AMP_ArrayVector
#define included_AMP_ArrayVector

#include <string>

#include "vectors/Vector.h"
#include "utils/FunctionTable.h"
#include "utils/Array.h"


namespace AMP {

namespace LinearAlgebra {

/** \brief A core-local vector
  * \details This is a Vector that implements the Vector interface for a std::vector<double>.
  */
template <typename T, typename FUN = FunctionTable, typename Allocator = std::allocator<T>>
class ArrayVector : public Vector
{
private:
    AMP::Array<T, FUN, Allocator> d_array;
    AMP_MPI d_comm;
    size_t d_globalSize;

    ArrayVector();
    ArrayVector( const ArrayVector & );

public:
    
    /** \brief    Create a ArrayVector
      * \details  This is the factory method for the ArrayVector.  It returns the shared pointer
      * to be used in the code
      * \param    localSize  The number of elements in the vector on this processor
      * \param    var The variable associated with the new vector
      */
    static Vector::shared_ptr create( const std::vector<size_t> &localSize,
                                      Variable::shared_ptr var );

    /** \brief    Create a ArrayVector
      * \details  This is the factory method for the ArrayVector.  It returns the shared pointer
      * to be used in the code
      * \param    localSize  The number of elements in the vector on this processor
      * \param    var The variable associated with the new vector
      * \param    comm The variable associated with the new vector
      */
    static Vector::shared_ptr
    create( const std::vector<size_t> &localSize, Variable::shared_ptr var, AMP_MPI comm );

    /** \brief    Create a ArrayVector
      * \details  This is the factory method for the ArrayVector.  It returns the shared pointer
      * to be used in the code that spans a comm and contains ghost values.
      * \param    var The variable associated with the new vector
      * \param    DOFs The DOFManager
      * \param    commlist The communication list
      */
    static Vector::shared_ptr create( Variable::shared_ptr var,
                                      AMP::Discretization::DOFManager::shared_ptr DOFs,
                                      AMP::LinearAlgebra::CommunicationList::shared_ptr commlist );

    /** \brief  Destructor
      */
    virtual ~ArrayVector() {}

    std::string type() const override { return "ArrayVector"; }


    /** \brief  Return the communicator this Vector spans
      */
    AMP_MPI getComm() const override { return d_comm; }

    using Vector::cloneVector;
    virtual Vector::shared_ptr cloneVector( const Variable::shared_ptr name ) const override;
    using Vector::copyVector;
    virtual void copyVector( Vector::const_shared_ptr src_vec ) override;
    virtual void swapVectors( Vector &other ) override;
    virtual void aliasVector( Vector &other ) override;

    //! resize the ArrayVector and reset the internal data structures
    void resize( const std::vector<size_t> &localDims );

    //! return a non-const reference to the internal data container
    Array<T, FUN, Allocator> &getArray( void ) { return d_array; }

    //! return a const reference to the internal data container
    const Array<T, FUN, Allocator> &getArray( void ) const { return d_array; }

    /** \brief Number of blocks of contiguous data in the Vector
      * \return Number of blocks in the Vector
      * \details  A vector is not necessarily contiguous in memory.  This method
      * returns the number of contiguous blocks in memory used by this vector
      */
    size_t numberOfDataBlocks() const override{ return 1; }

    /** \brief Number of elements in a data block
      * \param[in] i  particular data block
      * \return The size of a particular block
      */
    size_t sizeOfDataBlock( size_t i = 0 ) const override { NULL_USE(i); return d_array.length(); }

    /**\brief Copy data into this vector
      *\param[in] buf  Buffer to copy from
      */
    void putRawData( const double *buf ) override;

    /**\brief Copy data out of this vector
      *\param[out] buf  Buffer to copy to
      *\details The Vector should be pre-allocated to the correct size (getLocalSize())
      */
    void copyOutRawData( double *buf ) const override;

    /**\brief Number of elements "owned" by this core
      *\return  Number of entries stored contiguously on this processor
      *\details  For some types of variables, vectors may store "ghost"
      * data---possibly non-contiguous subsets of entries stored on other
      * cores.
      */
    size_t getLocalSize() const override { return d_array.length(); }

    /**\brief Number of total entries in this vector across all cores
      *\return Number of entries stored across all cores in this
      */
    size_t getGlobalSize() const override { return d_globalSize; }

    /**
      * \brief Set values in the vector by their local offset
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[in] vals the values to place in the vector
      * \details This will set the owned values for this core.  All indices are
      * from 0.
      * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
      */
    void setValuesByLocalID( int num, size_t *indices, const double *vals ) override { 
        NULL_USE(num);
        NULL_USE(indices);
        NULL_USE(vals);
        AMP_ERROR("Not implemented"); 
    }

    /**
      * \brief Set owned values using global identifier
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[in] vals the values to place in the vector
      *
      * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
      */
    void setLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) override { 
        NULL_USE(num);
        NULL_USE(indices);
        NULL_USE(vals);
        AMP_ERROR("Not implemented"); 
    }

    /**
      * \brief Add values to vector entities by their local offset
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[in] vals the values to place in the vector
      * \details This will set the owned values for this core.  All indices are
      * from 0.
      * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{this}_{\mathit{indices}_i} +
     * \mathit{vals}_i \f$
      */
    void addValuesByLocalID( int num, size_t *indices, const double *vals ) override { 
        NULL_USE(num);
        NULL_USE(indices);
        NULL_USE(vals);
        AMP_ERROR("Not implemented"); 
    }

    /**
      * \brief Add owned values using global identifier
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[in] vals the values to place in the vector
      *
      * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{this}_{\mathit{indices}_i} +
     * \mathit{vals}_i \f$
      */
    void addLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) override { 
        NULL_USE(num);
        NULL_USE(indices);
        NULL_USE(vals);
        AMP_ERROR("Not implemented"); 
    }

    /**
      * \brief Get local values in the vector by their global offset
      * \param[in] num  number of values to set
      * \param[in] indices the indices of the values to set
      * \param[out] vals the values to place in the vector
      * \details This will get any value owned by this core.
      */
    void getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const override  { 
        NULL_USE(num);
        NULL_USE(indices);
        NULL_USE(vals);
        AMP_ERROR("Not implemented"); 
    }

    /**\brief  A unique id for the underlying data allocation
      *\details This is a unique id that is associated with the data
      *   data allocation.  Views of a vector should preserve the id of
      *   the original vector.  Vectors that are not allocated, or contain
      *   multiple vectors (such as Multivector) should return 0.
      *   Note: this id is not consistent across multiple processors.
      */
    uint64_t getDataID() const override { return reinterpret_cast<uint64_t>( d_array.data() ); }

    /**
      * \brief This method is used to implement the assemble interface
      * of PETSc.
      * \details  This method is empty except for instantiations of NativePetscVector
      */
    void assemble() override { AMP_ERROR("Not implemented"); }

    /**
     * The next set of functions are derived from the VectorOperations class and should eventually
     * go away as we move away from a single type
     */

    /**
     * \param  alpha a scalar double
     * \brief  Set all compenents of a vector to a scalar.
     * For Vectors, the components of <em>this</em> are set to \f$\alpha\f$.
     */
    void setToScalar( double alpha ) override;

    /**
     * \param  alpha  a scalar double
     * \param  x  a vector
     * \brief  Set vector equal to scaled input.
     * For Vectors, \f$\mathit{this}_i = \alpha x_i\f$.
     */
    void scale( double alpha, const VectorOperations &x ) override;

    /**
     * \param  alpha  a scalar double
     *
     * \brief  Scale a vector.
     * For Vectors, \f$\mathit{this}_i = \alpha\mathit{this}_i\f$.
     */
    void scale( double alpha ) override;

    /**
     * \param  x  a vector
     * \param  y  a vector
     * \brief  Adds two vectors.
     * For Vectors, \f$\mathit{this}_i = x_i + y_i\f$.
     */
    void add( const VectorOperations &x, const VectorOperations &y ) override;

    /**
      * \param x  a vector
      * \param y  a vector
      * \brief Subtracts one vector from another.
      * For Vectors, \f$\mathit{this}_i = x_i - y_i\f$
     */
    void subtract( const VectorOperations &x, const VectorOperations &y ) override;

    /**
      * \param x  a vector
      * \param y  a vector
      * \brief Component-wise multiply one vector with another.
      * For Vectors, \f$\mathit{this}_i = x_i  y_i\f$
     */
    void multiply( const VectorOperations &x, const VectorOperations &y ) override;

    /**
      * \param x  a vector
      * \param y  a vector
      * \brief Component-wise divide one vector by another.
      * For Vectors, \f$\mathit{this}_i = x_i / y_i\f$
     */
    void divide( const VectorOperations &x, const VectorOperations &y ) override;

    /**
      * \param x  a vector
      * \brief Set this to the component-wise reciprocal of a vector.  \f$\mathit{this}_i =
     * 1/x_i\f$.
     */
    void reciprocal( const VectorOperations &x ) override;


    /**
     * \param alpha a scalar
     * \param x a vector
     * \param beta a scalar
     * \param y a vector
     * \brief Set a vector to be a linear combination of two vectors.
     *  \f$\mathit{this}_i = \alpha x_i + \beta y_i\f$.
     */
    void linearSum( double alpha,
                            const VectorOperations &x,
                            double beta,
                            const VectorOperations &y ) override;

    /**
      * \param alpha a scalar
      * \param x a vector
      * \param y a vector
      * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
     */
    void axpy( double alpha, const VectorOperations &x, const VectorOperations &y ) override;

    /**
      * \param alpha a scalar
      * \param beta a scalar
      * \param x  a vector
      * \brief Set this vector alpha * x + this.
      * \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
      */
    void axpby( double alpha, double beta, const VectorOperations &x ) override;

    /**
      * \param x a vector
      * \brief Set this to the component-wise absolute value of a vector.
      * \f$\mathit{this}_i = |x_i|\f$.
     */
    void abs( const VectorOperations &x ) override;

    /**
      * \brief Return the minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
     */
    double min( void ) const override;

    /**
      * \brief Return the maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
     */
    double max( void ) const override;

    /**
     * \brief Set data in this vector to random values on [0,1).
     */
    void setRandomValues( void ) override;
    using Vector::setRandomValues;

    /**
     * \brief Return discrete @f$ L_1 @f$ -norm of this vector.
     * \details Returns \f[\sum_i |\mathit{this}_i|\f]
     */
    double L1Norm( void ) const override;

    /**
     * \brief Return discrete @f$ L_2 @f$ -norm of this vector.
     * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
     */
    double L2Norm( void ) const override;

    /**
     * \brief Return the @f$ L_\infty @f$ -norm of this vector.
     * \details Returns \f[\max_i |\mathit{this}_i|\f]
     */
    double maxNorm( void ) const override;

    /**
      * \param x a vector
      * \brief Return the dot product of this vector with the argument vector.
      * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
     */
    using Vector::dot;
    double dot( const VectorOperations &x ) const override;

    /**
      * \fn equals (Vector & const rhs, double tol )
      * \brief  Determine if two vectors are equal using an absolute tolerance
      * \param[in] rhs Vector to compare to
      * \param[in] tol Tolerance of comparison
      * \return  True iff \f$||\mathit{rhs} - x||_\infty < \mathit{tol}\f$
      */
    bool equals( Vector const &rhs, double tol = 0.000001 ) const override;

protected:

    /** \brief Return a pointer to a particular block of memory in the
      * vector
      * \param i The block to return
      */
    void *getRawDataBlockAsVoid( size_t i ) override { AMP_ASSERT(i==0); return d_array.data(); }

    /** \brief Return a pointer to a particular block of memory in the
      * vector
      * \param i The block to return
      */
    const void *getRawDataBlockAsVoid( size_t i ) const override  { AMP_ASSERT(i==0); return d_array.data(); }
};

}
}

#include "ArrayVector.hpp"

#endif
