#ifndef included_AMP_MultiVectorData
#define included_AMP_MultiVectorData

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/vectors/data/DataChangeListener.h"
#include "AMP/vectors/data/VectorData.h"


namespace AMP::LinearAlgebra {


/**
  \brief  A class used to hold vector data

  \details

  MultiVectorData is a default implimentation of VectorData that stores
  the local values as a single block of data on the CPU.

  */
class MultiVectorData : public VectorData, public DataChangeListener
{

public: // Virtual functions
    //! Virtual destructor
    virtual ~MultiVectorData() {}


    //! Get the type name
    std::string VectorDataName() const override { return "MultiVectorData"; }


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
    void putRawData( const void *buf, const typeID &id ) override;

    /**\brief Copy data out of this vector
     *\param[out] buf  Buffer to copy to
     *\details The Vector should be pre-allocated to the correct size (getLocalSize())
     */
    void getRawData( void *buf, const typeID &id ) const override;


    void setValuesByLocalID( size_t, const size_t *, const void *, const typeID & ) override;
    void addValuesByLocalID( size_t, const size_t *, const void *, const typeID & ) override;
    void getValuesByLocalID( size_t, const size_t *, void *, const typeID & ) const override;
    void setGhostValuesByGlobalID( size_t, const size_t *, const void *, const typeID & ) override;
    void addGhostValuesByGlobalID( size_t, const size_t *, const void *, const typeID & ) override;
    void getGhostValuesByGlobalID( size_t, const size_t *, void *, const typeID & ) const override;
    size_t getGhostSize() const override;
    using VectorData::makeConsistent;
    void makeConsistent( ScatterType t ) override;
    UpdateState getUpdateStatus() const override;
    void setUpdateStatus( UpdateState state ) override;
    size_t getNumberOfComponents() const override;
    //! the next routine could be refined to depend on number of components
    bool hasContiguousData() const override { return false; }

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
    bool isType( const typeID &id, size_t block ) const override;

    /** \brief Swap the data with another VectorData object
     * \param rhs      The VectorData to swap with
     */
    void swapData( VectorData &rhs ) override;

    /** \brief Clone the data
     */
    std::shared_ptr<VectorData> cloneData() const override;

    void
    dumpOwnedData( std::ostream &out, size_t GIDoffset = 0, size_t LIDoffset = 0 ) const override;
    void dumpGhostedData( std::ostream &out, size_t offset = 0 ) const override;

    VectorData *getVectorData( size_t i );
    const VectorData *getVectorData( size_t i ) const;

    size_t getVectorDataSize() const { return d_data.size(); }

    AMP_MPI getComm() const override { return d_comm; }
    bool hasComm() const override { return true; }

    void assemble() override;

public:
    void receiveDataChanged() override { fireDataChange(); }

    explicit MultiVectorData( const AMP::AMP_MPI &comm )
        : d_comm( comm ), d_globalDOFManager( nullptr )
    {
    }

    void resetMultiVectorData( AMP::Discretization::DOFManager *manager,
                               const std::vector<VectorData *> &data );

protected:
    // Internal data
    AMP::AMP_MPI d_comm;
    std::vector<VectorData *> d_data;
    AMP::Discretization::DOFManager *d_globalDOFManager;
    std::vector<AMP::Discretization::DOFManager *> d_subDOFManager;


protected:
    /** A method that will translate an array of global ids relative to the multivector
     * into an array of arrays of global ids relative to the component vectors
     * \param[in] num            The number of DOFs that need to be mapped
     * \param[in] indices        The indices of the values relative to the multivector
     * \param[in] vals           Values associated somehow with the indices
     * \param[in] bytes          Size of a value
     * \param[out] out_indices   An array of arrays of mapped indices relative to constituent
     * vectors
     * \param[out] out_vals      The values partitioned according to out_indices
     * \param[out] remap         If not null, this is a list of where in the indices array an entry
     * comes from
     */
    void partitionGlobalValues( const int num,
                                const size_t *indices,
                                const void *vals,
                                const size_t bytes,
                                std::vector<std::vector<size_t>> &out_indices,
                                std::vector<std::vector<std::byte>> &out_vals,
                                std::vector<std::vector<int>> *remap = nullptr ) const;

    /** A method that will translate an array of local ids relative to the multivector
     * into an array of arrays of local ids relative to the component vectors
     * \param[in] num            The number of DOFs that need to be mapped
     * \param[in] indices        The indices of the values relative to the multivector
     * \param[in] vals           Values associated somehow with the indices
     * \param[in] bytes          Size of a value
     * \param[out] out_indices   An array of arrays of mapped indices relative to constituent
     * vectors
     * \param[out] out_vals      The values partitioned according to out_indices
     * \param[out] remap         If not null, this is a list of where in the indices array an entry
     * comes from
     */
    void partitionLocalValues( const int num,
                               const size_t *indices,
                               const void *vals,
                               const size_t bytes,
                               std::vector<std::vector<size_t>> &out_indices,
                               std::vector<std::vector<std::byte>> &out_vals,
                               std::vector<std::vector<int>> *remap = nullptr ) const;
};


} // namespace AMP::LinearAlgebra


#endif
