#ifndef included_AMP_VectorDataNull
#define included_AMP_VectorDataNull

#include "AMP/vectors/data/VectorData.h"


namespace AMP {
namespace LinearAlgebra {


template<typename TYPE>
class VectorDataIterator;


/**
  \brief  A class used to hold vector data

  \details

  VectorDataNull is a default implimentation of VectorData that stores
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
    inline void putRawData( const double * ) override {}

    /**\brief Copy data out of this vector
     *\details The Vector should be pre-allocated to the correct size (getLocalSize())
     */
    inline void copyOutRawData( double * ) const override {}

    /**\brief Number of elements "owned" by this core
      *\return  Number of entries stored contiguously on this processor
      *\details  For some types of variables, vectors may store "ghost"
      * data---possibly non-contiguous subsets of entries stored on other
      * cores.make

      */
    inline size_t getLocalSize() const override { return 0; }

    /**\brief Number of total entries in this vector across all cores
     *\return Number of entries stored across all cores in this
     */
    inline size_t getGlobalSize() const override { return 0; }

    /**
     * \brief Set values in the vector by their local offset
     * \param[in] num  number of values to set
     * \details This will set the owned values for this core.  All indices are
     * from 0.
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
     */
    inline void setValuesByLocalID( int num, size_t *, const double * ) override
    {
        AMP_INSIST( num == 0, "Cannot set values in NullVectorData" );
    }

    /**
     * \brief Set owned values using global identifier
     * \param[in] num  number of values to set
     *
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
     */
    inline void setLocalValuesByGlobalID( int num, size_t *, const double * ) override
    {
        AMP_INSIST( num == 0, "Cannot set values in NullVectorData" );
    }

    /**
     * \brief Add values to vector entities by their local offset
     * \param[in] num  number of values to set
     * \details This will set the owned values for this core.  All indices are
     * from 0.
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{this}_{\mathit{indices}_i} +
     * \mathit{vals}_i \f$
     */
    inline void addValuesByLocalID( int num, size_t *, const double * ) override
    {
        AMP_INSIST( num == 0, "Cannot add values in NullVectorData" );
    }

    /**
     * \brief Add owned values using global identifier
     * \param[in] num  number of values to set
     *
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{this}_{\mathit{indices}_i} +
     * \mathit{vals}_i \f$
     */
    inline void addLocalValuesByGlobalID( int num, size_t *, const double * ) override
    {
        AMP_INSIST( num == 0, "Cannot add values in NullVectorData" );
    }

    /**
     * \brief Get local values in the vector by their global offset
     * \param[in] num  number of values to set
     * \details This will get any value owned by this core.
     */
    inline void getLocalValuesByGlobalID( int num, size_t *, double * ) const override
    {
        AMP_INSIST( num == 0, "Cannot get values in NullVectorData" );
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
    inline bool isTypeId( size_t hash, size_t ) const override
    {
        return hash == typeid( TYPE ).hash_code();
    }

    inline void swapData( VectorData & ) override { AMP_ERROR( "Not finished" ); }

    inline std::shared_ptr<VectorData> cloneData() const override
    {
        return std::make_shared<VectorDataNull>();
    }
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
