#ifndef included_AMP_EpetraVectorData
#define included_AMP_EpetraVectorData

#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include "AMP/vectors/data/VectorData.h"

#include <Epetra_Vector.h>


namespace AMP {
namespace LinearAlgebra {

/** \class EpetraVectorEngineParameters
 * \brief Class that details how to construct an EpetraVectorEngine
 */
class EpetraVectorEngineParameters : public VectorParameters
{
public:
    /** \brief Constructor
     * \param[in] commList: communication list
     * \param[in] dofManager: DOFManager
     */
    EpetraVectorEngineParameters( std::shared_ptr<CommunicationList> commList,
                                  std::shared_ptr<AMP::Discretization::DOFManager> dofManager );
    /** \brief Constructor
        \param[in] local_size     The number of elements on this core
        \param[in] global_size    The number of elements in total
        \param[in] comm         Communicator to create the vector on
        \details  This assumes a contiguous allocation of data.  Core 0 has global ids
       \f$(0,1,\ldots,n-1)\f$, core 1
       has global ids \f$(n,n+1,n+2,\ldots,m)\f$, etc.
        */
    EpetraVectorEngineParameters( size_t local_size, size_t global_size, const AMP_MPI &comm );

    /** \brief Constructor
     * \param[in]  local_size    The number of elements on this core
     * \param[in]  global_size   The number of elements in total
     * \param[in]  emap        An Epetra_Map for the data
     * \param[in]  ecomm       An Epetra_MpiComm for constructing the vector on
     * \details  This allows construction of an EpetraVectorEngine from handy Epetra objects
     */
    EpetraVectorEngineParameters( size_t local_size,
                                  size_t global_size,
                                  std::shared_ptr<Epetra_Map> emap,
                                  const AMP_MPI &ecomm );

    //! Destructor
    virtual ~EpetraVectorEngineParameters();

    /** \brief  Return the Epetra_Map for this engine
     * \return  The Epetra_Map
     */
    Epetra_Map &getEpetraMap();

    //! Return the local size
    inline size_t getLocalSize() const { return d_end - d_begin; }

    //! Return the local size
    inline size_t getGlobalSize() const { return d_global; }

    //! Return the first DOF on this core
    inline size_t beginDOF() const { return d_begin; }

    /** \brief  Return the Epetra_MpiComm for this engine
     * \return The Epetra_MpiComm
     */
    inline AMP_MPI getComm() const { return d_comm; }

private:
    size_t d_begin  = 0;                          // Starting DOF
    size_t d_end    = 0;                          // Ending DOF
    size_t d_global = 0;                          // Number of global DOFs
    AMP_MPI d_comm;                               // Comm
    std::shared_ptr<Epetra_Map> d_emap = nullptr; // Epetra map
};

/** \class EpetraVectorData
 * \brief A linear algebra engine that uses Epetra
 * \details  Use the Epetra implementation of the L1 BLAS routines.  Unlike other
 * libraries, it is very difficult to separate the data from the engine.  For this
 * reason, the EpetraVectorEngine contains the Epetra_Vector to operate on.
 */
class EpetraVectorData : public VectorData
{

public: // Virtual functions
    //! Virtual destructor
    virtual ~EpetraVectorData() {}
    static EpetraVectorData *create( std::shared_ptr<EpetraVectorEngineParameters> alias,
				     std::shared_ptr<VectorData> buf );

    std::string VectorDataName() const override { return "EpetraVectorData"; }
    size_t getLocalSize() const override { return d_iLocalSize; }
    size_t getGlobalSize() const override { return d_iGlobalSize; }
    size_t getLocalStartID() const override { return d_iLocalStart; }
    size_t numberOfDataBlocks() const override { return 1; }
    size_t sizeOfDataBlock( size_t i ) const override { return i == 0 ? d_iLocalSize : 0; }
    void setValuesByLocalID( int i, size_t *, const double *val ) override;
    void setLocalValuesByGlobalID( int i, size_t *, const double *val ) override;
    void addValuesByLocalID( int i, size_t *, const double *val ) override;
    void addLocalValuesByGlobalID( int i, size_t *, const double *val ) override;
    void getValuesByLocalID( int i, size_t *, double *val ) const override;
    void getLocalValuesByGlobalID( int i, size_t *, double *val ) const override;
    void putRawData( const double *in ) override;
    void copyOutRawData( double *out ) const override;
    uint64_t getDataID() const override
    {
        return reinterpret_cast<uint64_t>( getRawDataBlockAsVoid( 0 ) );
    }
    void *getRawDataBlockAsVoid( size_t i ) override;
    const void *getRawDataBlockAsVoid( size_t i ) const override;
    size_t sizeofDataBlockType( size_t ) const override { return sizeof( double ); }
    bool isTypeId( size_t hash, size_t ) const override
    {
        return hash == typeid( double ).hash_code();
    }
    void swapData( VectorData & ) override { AMP_ERROR( "Not finished" ); }
    std::shared_ptr<VectorData> cloneData() const override;

    /** \brief  Get the raw Epetra_Vector
     * \return  The Epetra_Vector currently used by this engine
     */
    inline Epetra_Vector &getEpetra_Vector() { return d_epetraVector; }

    /** \brief  Get the raw Epetra_Vector
     * \return  The Epetra_Vector currently used by this engine
     */
    inline const Epetra_Vector &getEpetra_Vector() const { return d_epetraVector; }
 
 protected:
   EpetraVectorData( std::shared_ptr<EpetraVectorEngineParameters> alias,
		     Epetra_DataAccess,
		     const Epetra_BlockMap &,
		     std::shared_ptr<VectorData>,
		     int, int, int );

    //! The Epetra_Vector to perform work on
    Epetra_Vector d_epetraVector;

    //! A shared ptr to the buffer (to prevent it from leaving scope)
    std::shared_ptr<VectorData> d_buf_scope;

    //! The local start index
    int d_iLocalStart;

    //! The number of local elements in the vector
    int d_iLocalSize;

    //! The number of elements in the entire vector
    int d_iGlobalSize;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
