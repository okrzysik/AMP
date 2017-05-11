#ifndef included_AMP_ArrayVector
#define included_AMP_ArrayVector

#include <string>

#include "vectors/Vector.h"
#include "vectors/operations/VectorOperationsDefault.h"
#include "vectors/operations/VectorOperationsDefault.hpp"
#include "utils/FunctionTable.h"
#include "utils/Array.h"


namespace AMP {

namespace LinearAlgebra {

/** \brief A core-local vector
  * \details This is a Vector that implements the Vector interface for a std::vector<double>.
  */
template <typename T, typename FUN = FunctionTable, typename Allocator = std::allocator<T>>
class ArrayVector :
    public Vector,
    public VectorOperationsDefault<T>
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

    virtual bool isTypeId( size_t hash, size_t ) const override { return hash == typeid(T).hash_code(); }

};

}
}

#include "ArrayVector.hpp"

#endif
