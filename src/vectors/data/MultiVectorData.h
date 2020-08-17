#ifndef included_AMP_MultiVectorData
#define included_AMP_MultiVectorData

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/vectors/DataChangeListener.h"
#include "AMP/vectors/data/VectorData.h"


namespace AMP {
namespace LinearAlgebra {


/**
  \brief  A class used to hold vector data

  \details

  MultiVectorData is a default implimentation of VectorData that stores
  the local values as a single block of data on the CPU.

  */
class MultiVectorData : virtual public VectorData, public DataChangeListener
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
      * cores.make

      */
    size_t getLocalSize() const override;

    /**\brief Number of total entries in this vector across all cores
     *\return Number of entries stored across all cores in this
     */
    size_t getGlobalSize() const override;

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

    void setGhostValuesByGlobalID( int num, size_t *indices, const double *in_vals ) override;
    void setValuesByGlobalID( int num, size_t *indices, const double *in_vals ) override;
    void addValuesByGlobalID( int num, size_t *indices, const double *in_vals ) override;
    void getValuesByGlobalID( int num, size_t *indices, double *out_vals ) const override;
    void getGhostValuesByGlobalID( int num, size_t *indices, double *out_vals ) const override;
    void getValuesByLocalID( int num, size_t *indices, double *out_vals ) const override;
    size_t getGhostSize() const override;
    void makeConsistent( ScatterType t ) override;
    UpdateState getUpdateStatus() const override;
    void setUpdateStatus( UpdateState state ) override;


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

    VectorData *getVectorData( size_t i );
    const VectorData *getVectorData( size_t i ) const;

    size_t numberOfComponents( void ) const { return d_data.size(); }


public:
    void recieveDataChanged() override { fireDataChange(); }


protected:
    MultiVectorData() : d_globalDOFManager( nullptr ) {}

    // Internal data
    std::vector<VectorData *> d_data;
    AMP::Discretization::DOFManager *d_globalDOFManager;
    std::vector<AMP::Discretization::DOFManager *> d_subDOFManager;


protected:
    /** A method that will translate an array of global ids relative to the multivector
     * into an array of arrays of global ids relative to the component vectors
     * \param[in] num            The number of DOFs that need to be mapped
     * \param[in] indices        The indices of the values relative to the multivector
     * \param[in] vals           Values associated somehow with the indices
     * \param[out] out_indices   An array of arrays of mapped indices relative to constituent
     * vectors
     * \param[out] out_vals      The values partitioned according to out_indices
     * \param[out] remap         If not null, this is a list of where in the indices array an entry
     * comes from
     */
    void partitionGlobalValues( const int num,
                                const size_t *indices,
                                const double *vals,
                                std::vector<std::vector<size_t>> &out_indices,
                                std::vector<std::vector<double>> &out_vals,
                                std::vector<std::vector<int>> *remap = nullptr ) const;

    /** A method that will translate an array of local ids relative to the multivector
     * into an array of arrays of local ids relative to the component vectors
     * \param[in] num            The number of DOFs that need to be mapped
     * \param[in] indices        The indices of the values relative to the multivector
     * \param[in] vals           Values associated somehow with the indices
     * \param[out] out_indices   An array of arrays of mapped indices relative to constituent
     * vectors
     * \param[out] out_vals      The values partitioned according to out_indices
     * \param[out] remap         If not null, this is a list of where in the indices array an entry
     * comes from
     */
    void partitionLocalValues( const int num,
                               const size_t *indices,
                               const double *vals,
                               std::vector<std::vector<size_t>> &out_indices,
                               std::vector<std::vector<double>> &out_vals,
                               std::vector<std::vector<int>> *remap = nullptr ) const;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
