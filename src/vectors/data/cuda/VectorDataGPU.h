#ifndef included_AMP_VectorDataGPU
#define included_AMP_VectorDataGPU

#include "AMP/vectors/data/VectorData.h"


namespace AMP {
namespace LinearAlgebra {


template<typename TYPE>
class VectorDataIterator;


/**
  \brief  A class used to hold vector data

  \details

  VectorDataGPU is a default implimentation of VectorData that stores
  the local values as a single block of data on the GPU.

  */
template<typename TYPE = double>
class VectorDataGPU : public VectorData
{
public: // Constructors
    VectorDataGPU( size_t start, size_t localSize, size_t globalSize );


public: // Virtual functions
    //! Virtual destructor
    virtual ~VectorDataGPU();

    //! Get the type name
    virtual std::string VectorDataName() const override;

    /** \brief Number of blocks of contiguous data in the Vector
     * \return Number of blocks in the Vector
     * \details  A vector is not necessarily contiguous in memory.  This method
     * returns the number of contiguous blocks in memory used by this vector
     */
    virtual size_t numberOfDataBlocks() const override;

    /** \brief Number of elements in a data block
     * \param[in] i  particular data block
     * \return The size of a particular block
     */
    virtual size_t sizeOfDataBlock( size_t i = 0 ) const override;


    /**\brief Copy data into this vector
     *\param[in] buf  Buffer to copy from
     */
    virtual void putRawData( const double *buf ) override;

    /**\brief Copy data out of this vector
     *\param[out] buf  Buffer to copy to
     *\details The Vector should be pre-allocated to the correct size (getLocalSize())
     */
    virtual void copyOutRawData( double *buf ) const override;

    /**\brief Number of elements "owned" by this core
      *\return  Number of entries stored contiguously on this processor
      *\details  For some types of variables, vectors may store "ghost"
      * data---possibly non-contiguous subsets of entries stored on other
      * cores.make

      */
    virtual size_t getLocalSize() const override;

    /**\brief Number of total entries in this vector across all cores
     *\return Number of entries stored across all cores in this
     */
    virtual size_t getGlobalSize() const override;

    /**
     * \brief Set values in the vector by their local offset
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     * \details This will set the owned values for this core.  All indices are
     * from 0.
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
     */
    virtual void setValuesByLocalID( int num, size_t *indices, const double *vals ) override;

    /**
     * \brief Set owned values using global identifier
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     *
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
     */
    virtual void setLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) override;

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
    virtual void addValuesByLocalID( int num, size_t *indices, const double *vals ) override;

    /**
     * \brief Add owned values using global identifier
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     *
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{this}_{\mathit{indices}_i} +
     * \mathit{vals}_i \f$
     */
    virtual void addLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) override;

    /**
     * \brief Get local values in the vector by their global offset
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[out] vals the values to place in the vector
     * \details This will get any value owned by this core.
     */
    virtual void getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const override;


public: // Advanced virtual functions
    /**\brief  A unique id for the underlying data allocation
     *\details This is a unique id that is associated with the data
     *   data allocation.  Views of a vector should preserve the id of
     *   the original vector.  Vectors that are not allocated, or contain
     *   multiple vectors (such as Multivector) should return 0.
     *   Note: this id is not consistent across multiple processors.
     */
    virtual uint64_t getDataID() const override;

    /** \brief Return a pointer to a particular block of memory in the vector
     * \param i The block to return
     */
    virtual void *getRawDataBlockAsVoid( size_t i ) override;

    /** \brief Return a pointer to a particular block of memory in the
     * vector
     * \param i        The block to return
     */
    virtual const void *getRawDataBlockAsVoid( size_t i ) const override;

    /** \brief Return the result of sizeof(TYPE) for the given data block
     * \param i The block to return
     */
    virtual size_t sizeofDataBlockType( size_t i ) const override;

    /** \brief Is the data of the given type
     * \param hash     The hash code: typeid(myint).hash_code()
     * \param block    The block id to check
     */
    virtual bool isTypeId( size_t hash, size_t block ) const override;

    /** \brief Swap the data with another VectorData object
     * \param rhs      The VectorData to swap with
     */
    virtual void swapData( VectorData &rhs ) override;

    /** \brief Clone the data
     */
    virtual std::shared_ptr<VectorData> cloneData() const override;


public: // Non-virtual functions
    /** \brief Access the raw element
     * \param i        The element to return (local index)
     */
    TYPE &operator[]( size_t i );

    /** \brief Access the raw element
     * \param i        The element to return (local index)
     */
    const TYPE &operator[]( size_t i ) const;


protected:
    VectorDataGPU() : d_Data( nullptr ), d_startIndex( 0 ), d_localSize( 0 ), d_globalSize( 0 ) {}

    void allocate( size_t start, size_t localSize, size_t globalSize );


private:
    TYPE *d_Data;
    size_t d_startIndex;
    size_t d_localSize;
    size_t d_globalSize;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
