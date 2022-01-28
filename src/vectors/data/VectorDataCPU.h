#ifndef included_AMP_VectorDataCPU
#define included_AMP_VectorDataCPU

#include "AMP/utils/UtilityMacros.h"
#include "AMP/vectors/data/VectorData.h"


namespace AMP::LinearAlgebra {


template<typename TYPE>
class VectorDataIterator;


/**
  \brief  A class used to hold vector data

  \details

  VectorDataCPU is a default implimentation of VectorData that stores
  the local values as a single block of data on the CPU.

  */
template<typename TYPE = double>
class VectorDataCPU : public VectorData
{
public: // Constructors
    VectorDataCPU( size_t start, size_t localSize, size_t globalSize );

    VectorDataCPU( const VectorDataCPU & ) = delete;

public: // Virtual functions
    //! Virtual destructor
    virtual ~VectorDataCPU() {}

    //! Get the type name
    std::string VectorDataName() const override;

    /** \brief Number of blocks of contiguous data in the Vector
     * \return Number of blocks in the Vector
     * \details  A vector is not necessarily contiguous in memory.  This method
     * returns the number of contiguous blocks in memory used by this vector
     */
    size_t numberOfDataBlocks() const override;

    /** \brief Number of elements in a data block
     * \param[in] i  particular data block
     * \return The size of a particular block
     */
    size_t sizeOfDataBlock( size_t i = 0 ) const override;


    /**\brief Copy data into this vector
     *\param[in] buf  Buffer to copy from
     */
    void putRawData( const double *buf ) override;

    /**\brief Copy data out of this vector
     *\param[out] buf  Buffer to copy to
     *\details The Vector should be pre-allocated to the correct size (getLocalSize())
     */
    void copyOutRawData( double *buf ) const override;

    /**
     * \brief Set values in the vector by their local offset
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     * \details This will set the owned values for this core.  All indices are
     * from 0.
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
     */
    void setValuesByLocalID( int num, size_t *indices, const double *vals ) override;

    /**
     * \brief Set owned values using global identifier
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     *
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
     */
    void setLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) override;

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
    void addValuesByLocalID( int num, size_t *indices, const double *vals ) override;

    /**
     * \brief Add owned values using global identifier
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     *
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{this}_{\mathit{indices}_i} +
     * \mathit{vals}_i \f$
     */
    void addLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) override;

    /**
     * \brief Get local values in the vector by their global offset
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[out] vals the values to place in the vector
     * \details This will get any value owned by this core.
     */
    void getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const override;


public: // Advanced virtual functions
    /**\brief  A unique id for the underlying data allocation
     *\details This is a unique id that is associated with the data
     *   data allocation.  Views of a vector should preserve the id of
     *   the original vector.  Vectors that are not allocated, or contain
     *   multiple vectors (such as Multivector) should return 0.
     *   Note: this id is not consistent across multiple processors.
     */
    uint64_t getDataID() const override;

    /** \brief Return a pointer to a particular block of memory in the vector
     * \param i The block to return
     */
    void *getRawDataBlockAsVoid( size_t i ) override;

    /** \brief Return a pointer to a particular block of memory in the
     * vector
     * \param i        The block to return
     */
    const void *getRawDataBlockAsVoid( size_t i ) const override;

    /** \brief Return the result of sizeof(TYPE) for the given data block
     * \param i The block to return
     */
    size_t sizeofDataBlockType( size_t i ) const override;

    /** \brief Is the data of the given type
     * \param hash     The hash code: typeid(myint).hash_code()
     * \param block    The block id to check
     */
    bool isTypeId( size_t hash, size_t block ) const override;

    /** \brief Swap the data with another VectorData object
     * \param rhs      The VectorData to swap with
     */
    void swapData( VectorData &rhs ) override;

    /** \brief Clone the data
     */
    std::shared_ptr<VectorData> cloneData() const override;


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
    VectorDataCPU() {}

    void allocate( size_t start, size_t localSize, size_t globalSize );


private:
    std::vector<TYPE> d_Data;
};


} // namespace AMP::LinearAlgebra


#endif
