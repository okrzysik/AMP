#ifndef included_AMP_VectorData
#define included_AMP_VectorData

#include "AMP/vectors/CommunicationList.h"
#include "AMP/vectors/data/DataChangeFirer.h"

#include <memory>
#include <vector>


namespace AMP {
namespace LinearAlgebra {


template<typename TYPE>
class VectorDataIterator;


/**
 * \brief  A class used to hold vector data
 * \details  VectorData is a class to helping disassociate data storage
 *    and vector operations such as dot product, norms, etc.
 */
class VectorData : public DataChangeFirer
{
public: // enums
    /**\brief Flag to choose algorithm for makeConsistent
     *\see makeConsistent
     */
    enum class ScatterType { CONSISTENT_ADD, CONSISTENT_SET };

    /**\brief The four states a Vector can be in
     *\see makeConsistent
     */
    enum class UpdateState { UNCHANGED, LOCAL_CHANGED, ADDING, SETTING, MIXED };


public: // Virtual functions
    //! Virtual destructor
    virtual ~VectorData() {}


    //! Get the type name
    virtual std::string VectorDataName() const = 0;


    /** \brief Number of blocks of contiguous data in the Vector
     * \return Number of blocks in the Vector
     * \details  A vector is not necessarily contiguous in memory.  This method
     * returns the number of contiguous blocks in memory used by this vector
     */
    virtual size_t numberOfDataBlocks() const = 0;

    /** \brief Number of elements in a data block
     * \param[in] i  particular data block
     * \return The size of a particular block
     */
    virtual size_t sizeOfDataBlock( size_t i = 0 ) const = 0;


    /**\brief Copy data into this vector
     *\param[in] buf  Buffer to copy from
     */
    virtual void putRawData( const double *buf ) = 0;

    /**\brief Copy data out of this vector
     *\param[out] buf  Buffer to copy to
     *\details The Vector should be pre-allocated to the correct size (getLocalSize())
     */
    virtual void copyOutRawData( double *buf ) const = 0;

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

    /**\brief get local start id core.
     *\return The first entry "owned" by this core
     */
    virtual size_t getLocalStartID() const;

    /**\brief Number of entries "owned" by other cores stored on this
     * core.
     *\return Number of entries "owned" by other cores stored on this core
     */
    virtual size_t getGhostSize() const;

    /**
     * \brief Set values in the vector by their local offset
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     * \details This will set the owned values for this core.  All indices are
     * from 0.
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
     */
    virtual void setValuesByLocalID( int num, size_t *indices, const double *vals ) = 0;


    /**
     * \brief Set owned values using global identifier
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     *
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
     */
    virtual void setLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) = 0;

    /**
     * \brief Set ghost values using global identifier
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     *
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{vals}_i \f$
     */
    virtual void setGhostValuesByGlobalID( int num, size_t *indices, const double *vals );

    /**
     * \brief Set owned or shared values using global identifier
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     * \details Since the shared buffer and owned buffer are separate,
     * this function must sort the data by buffer before setting
     * values.
     */
    virtual void setValuesByGlobalID( int num, size_t *indices, const double *vals );

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
    virtual void addValuesByLocalID( int num, size_t *indices, const double *vals ) = 0;

    /**
     * \brief Add owned values using global identifier
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     *
     * \f$ \mathit{this}_{\mathit{indices}_i} = \mathit{this}_{\mathit{indices}_i} +
     * \mathit{vals}_i \f$
     */
    virtual void addLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) = 0;

    /**
     * \brief Add owned or shared values using global identifier
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     * \details Since the shared buffer and owned buffer are separate,
     * this function must sort the data by buffer before setting
     * values.
     */
    virtual void addValuesByGlobalID( int num, size_t *indices, const double *vals );


    /**
     * \brief get ghosted values to add to off-proc elements
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[in] vals the values to place in the vector
     * \details This will get the ghosted updates this processor has made.  All indices are
     * from global 0.
     */
    virtual void getGhostAddValuesByGlobalID( int num, size_t *indices, double *vals ) const;

    /**
     * \brief get values in the vector by their local offset
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[out] vals the values to place in the vector
     * \details This will get the owned values for this core.  All indices are
     * from 0.
     */
    virtual void getValuesByGlobalID( int num, size_t *indices, double *vals ) const;

    /**
     * \brief Get local values in the vector by their global offset
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[out] vals the values to place in the vector
     * \details This will get any value owned by this core.
     */
    virtual void getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const = 0;

    /**
     * \brief Get ghost values in the vector by their global offset
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[out] vals the values to place in the vector
     * \details This will get any value owned by this core.
     */
    virtual void getGhostValuesByGlobalID( int num, size_t *indices, double *vals ) const;

    /**
     * \brief Get local values in the vector by their global offset
     * \param[in] num  number of values to set
     * \param[in] indices the indices of the values to set
     * \param[out] vals the values to place in the vector
     * \details This will get any value used by this core.
     */
    virtual void getValuesByLocalID( int num, size_t *indices, double *vals ) const;


public: // Advanced functions
    /**\brief  A unique id for the underlying data allocation
     *\details This is a unique id that is associated with the data
     *   data allocation.  Views of a vector should preserve the id of
     *   the original vector.  Vectors that are not allocated, or contain
     *   multiple vectors (such as Multivector) should return 0.
     *   Note: this id is not consistent across multiple processors.
     */
    virtual uint64_t getDataID() const = 0;

    /** \brief Return a pointer to a particular block of memory in the vector
     * \param i The block to return
     */
    virtual void *getRawDataBlockAsVoid( size_t i ) = 0;

    /** \brief Return a pointer to a particular block of memory in the
     * vector
     * \param i The block to return
     */
    virtual const void *getRawDataBlockAsVoid( size_t i ) const = 0;

    /** \brief Return the result of sizeof(TYPE) for the given data block
     * \param i The block to return
     */
    virtual size_t sizeofDataBlockType( size_t i ) const = 0;

    /** \brief Is the data of the given type
     * \param hash     The hash code: typeid(myint).hash_code()
     * \param block    The block id to check
     */
    virtual bool isTypeId( size_t hash, size_t block ) const = 0;

    /** \brief Swap the data with another VectorData object
     * \param rhs      The VectorData to swap with
     */
    virtual void swapData( VectorData &rhs ) = 0;

    /** \brief Clone the data
     */
    virtual std::shared_ptr<VectorData> cloneData() const = 0;

    /** \brief Associate the ghost buffer of a Vector with this Vector
     * \param in  The Vector to share a ghost buffer with
     */
    void aliasGhostBuffer( std::shared_ptr<VectorData> in );

    /** \brief Check if the two VectorData objects are alias of each other
     * \details  This function checks if two VectorData objects are alias of each other.
     *     Two VectorData objects are alias if their data blocks are the same size and
     *     point to the same memory blocks.
     * \param[in] rhs  VectorData to compare
     * \see makeConsistent
     */
    virtual bool isAnAliasOf( const VectorData &rhs ) const;

public:
    /** \brief Write owned data to an std::ostream
     * \param[in] out  The output stream to write to.
     * \param[in] GIDoffset  A number to add to the global ID when writing information
     * \param[in] LIDoffset  A number to add to the local ID when writing information
     */
    virtual void
    dumpOwnedData( std::ostream &out, size_t GIDoffset = 0, size_t LIDoffset = 0 ) const;

    /** \brief Write data owned by other processors to an std::ostream
     * \param[in] out  The output stream to write to.
     * \param[in] offset  A number to add to the global ID when writing information
     */
    virtual void dumpGhostedData( std::ostream &out, size_t offset = 0 ) const;


public: // Virtual functions dealing with the update status
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
    virtual UpdateState getUpdateStatus() const;


    /** \brief  Sets the current update state of the Vector
     * \details  This sets the update status of the vector.
     * This function should only be called by advanced users
     * \param[in] state  State of the vector to set
     */
    virtual void setUpdateStatus( UpdateState state );

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

    //! Get the communicator
    virtual AMP_MPI getComm() const;

    virtual bool hasComm( void ) const { return ( d_CommList != nullptr ); }

    //! Get the CommunicationList for this Vector
    virtual CommunicationList::shared_ptr getCommunicationList() const;

    /**\brief Set the CommunicationList for this Vector
     *\details  Setting the CommunicationList for a Vector may involve
     * reallocating ghost storage.
     */
    virtual void setCommunicationList( CommunicationList::shared_ptr comm );

public: // Non-virtual functions
    /**
     * \brief Return an iterator to the beginning of the data
     * \returns A VectorDataIterator
     * \details Since the Vector presents an interface to a contiguous
     *     block of data, it is natural for it to provide a random
     *     access iterator.
     * \warning The non-const version of the iterators will automatically
     *     leave the vector in a non-consistent state.  The user may
     *     be required to call makeConsistent.
     */
    template<class TYPE = double>
    inline VectorDataIterator<TYPE> begin();

    /// @copydoc VectorData::begin()
    template<class TYPE = double>
    inline VectorDataIterator<const TYPE> begin() const;

    /// @copydoc VectorData::begin()
    template<class TYPE = double>
    inline VectorDataIterator<const TYPE> constBegin() const;

    /**
     * \brief Return an iterator to the end of the data
     * \returns A VectorDataIterator
     * \details Since the Vector presents an interface to a contiguous
     *     block of data, it is natural for it to provide a random
     *     access iterator.
     * \warning The non-const version of the iterators will automatically
     *     leave the vector in a non-consistent state.  The user may
     *     be required to call makeConsistent.
     */
    template<class TYPE = double>
    inline VectorDataIterator<TYPE> end();

    /// @copydoc VectorData::end()
    template<class TYPE = double>
    inline VectorDataIterator<const TYPE> end() const;

    /// @copydoc VectorData::end()
    template<class TYPE = double>
    inline VectorDataIterator<const TYPE> constEnd() const;

    /** \brief Obtain a particular contiguous block of data cast to RETURN_TYPE
     * \tparam RETURN_TYPE  The pointer type of the return
     * \param[in] i  Which block
     * \return A contiguous array of type RETURN_TYPE
     */
    template<typename RETURN_TYPE>
    RETURN_TYPE *getRawDataBlock( size_t i = 0 );

    /** \brief Obtain a particular contiguous block of data cast to RETURN_TYPE
     * \tparam RETURN_TYPE  The pointer type of the return
     * \param[in] i  Which block
     * \return A const contiguous array of type RETURN_TYPE
     */
    template<typename RETURN_TYPE>
    const RETURN_TYPE *getRawDataBlock( size_t i = 0 ) const;

    /** \brief Check if the data is of the given type
     * \return true if all data is of the given type
     */
    template<typename TYPE>
    bool isType() const;

    /** \brief Check if the data is of the given type
     * \return true if the ith block of data is of the given type
     */
    template<typename TYPE>
    bool isBlockType( size_t i = 0 ) const;

    /** \brief Copy ghosted vlues to a vector
     * \param[in] rhs  Vector to copy ghost values from
     * \details  This ensures that ghosted values on this and rhs
     * are the same without a call to makeConsistent.
     * \see makeConsistent
     */
    void copyGhostValues( const VectorData &rhs );

public:
    /** \brief Notify listeners that data has changed in this vector.
     */
    virtual void dataChanged();

    /* \brief  Returns true if this vector has this element
     * \param[in]  GID  The global ID of the element
     */
    virtual bool containsGlobalElement( size_t GID );

    /**
     * \brief This method is used to implement the assemble interface
     * of PETSc.
     * \details  This method is empty except for instantiations of NativePetscVectorData
     */
    virtual void assemble() {}

public: // Non virtual functions
    /** \brief  Return the current update state of this Vector
     * \details  This returns the pointer to the update state
     *  of the current vector only (not vectors it contains).
     *  It should NOT be used by users.
     */
    std::shared_ptr<UpdateState> getUpdateStatusPtr() const;

    /** \brief  Tie the current update state to another
     * \details  This sets the pointer to the update state
     *  of the current vector only (not vectors it contains).
     *  It should NOT be used by users.
     * \param  rhs Pointer to share update state with
     */
    void setUpdateStatusPtr( std::shared_ptr<UpdateState> rhs );

    bool hasGhosts( void ) { return ( d_Ghosts != nullptr ); }

    std::vector<double> &getGhosts() { return *d_Ghosts; }

protected: // Internal data
    //! The communication list for this vector
    CommunicationList::shared_ptr d_CommList = nullptr;

    /** \brief  The current update state for a vector
     * \details A Vector can be in one of three states.
     *  This is the current state of the vector
     *  Because a vector can be composed of vectors,
     *  the update state needs to be shared between them.
     */
    std::shared_ptr<UpdateState> d_UpdateState = nullptr;

    // Ghost data
    std::shared_ptr<std::vector<double>> d_Ghosts    = nullptr;
    std::shared_ptr<std::vector<double>> d_AddBuffer = nullptr;

    // Friends
    friend class VectorOperations;


public:
    //! Default constructor
    VectorData();

    VectorData( std::shared_ptr<CommunicationList> commList );
};


} // namespace LinearAlgebra
} // namespace AMP


#include "AMP/vectors/data/VectorData.inline.h"


#endif
