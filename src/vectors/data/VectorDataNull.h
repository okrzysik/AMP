#ifndef included_AMP_VectorDataNull
#define included_AMP_VectorDataNull

#include "AMP/vectors/data/VectorData.h"


namespace AMP::LinearAlgebra {


template<typename TYPE>
class VectorDataIterator;


/**
  \brief  A class used to hold vector data

  \details

  VectorDataNull is a default implementation of VectorData that stores
  the local values as a single block of data on the CPU.

  */
template<typename TYPE = double>
class VectorDataNull : public VectorData
{
public:
    VectorDataNull() {}

public: // Virtual functions
    //! Virtual destructor
    virtual ~VectorDataNull() {}


    //! Get the type name
    virtual std::string VectorDataName() const override { return "VectorDataNull"; }


    /** \brief Number of blocks of contiguous data in the Vector
     * \return Number of blocks in the Vector
     * \details  A vector is not necessarily contiguous in memory.  This method
     * returns the number of contiguous blocks in memory used by this vector
     */
    inline size_t numberOfDataBlocks() const override { return 0; }

    /** \brief Number of elements in a data block
     * \return The size of a particular block
     */
    inline size_t sizeOfDataBlock( size_t = 0 ) const override { return 0; }


    /**\brief Copy data into this vector
     */
    inline void putRawData( const void *, const typeID & ) override {}

    /**\brief Copy data out of this vector
     *\details The Vector should be pre-allocated to the correct size (getLocalSize())
     */
    inline void getRawData( void *, const typeID & ) const override {}

    inline void
    getValuesByLocalID( size_t N, const size_t *, void *, const typeID & ) const override
    {
        AMP_INSIST( N == 0, "Cannot get values in NullVectorData" );
    }
    inline void
    setValuesByLocalID( size_t N, const size_t *, const void *, const typeID & ) override
    {
        AMP_INSIST( N == 0, "Cannot set values in NullVectorData" );
    }
    inline void
    addValuesByLocalID( size_t N, const size_t *, const void *, const typeID & ) override
    {
        AMP_INSIST( N == 0, "Cannot add values in NullVectorData" );
    }


public: // Advanced virtual functions
    /**\brief  A unique id for the underlying data allocation
     *\details This is a unique id that is associated with the data
     *   data allocation.  Views of a vector should preserve the id of
     *   the original vector.  Vectors that are not allocated, or contain
     *   multiple vectors (such as Multivector) should return 0.
     *   Note: this id is not consistent across multiple processors.
     */
    inline uint64_t getDataID() const override { return 0; }

    /** \brief Return a pointer to a particular block of memory in the vector
     */
    inline void *getRawDataBlockAsVoid( size_t ) override { return nullptr; }

    /** \brief Return a pointer to a particular block of memory in the
     * vector
     */
    inline const void *getRawDataBlockAsVoid( size_t ) const override { return nullptr; }

    /** \brief Return the result of sizeof(TYPE) for the given data block
     */
    inline size_t sizeofDataBlockType( size_t ) const override { return sizeof( TYPE ); }

    /** \brief Is the data of the given type
     * \param hash     The hash code: typeid(myint).hash_code()
     */
    inline bool isType( const typeID &id, size_t ) const override
    {
        return id == getTypeID<TYPE>();
    }

    inline void swapData( VectorData & ) override { AMP_ERROR( "Not finished" ); }

    inline std::shared_ptr<VectorData> cloneData() const override
    {
        return std::make_shared<VectorDataNull>();
    }
};


} // namespace AMP::LinearAlgebra


#endif
