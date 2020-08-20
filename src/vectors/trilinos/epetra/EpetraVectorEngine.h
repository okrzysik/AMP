#ifndef included_AMP_EpetraVectorEngine
#define included_AMP_EpetraVectorEngine

#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include "AMP/vectors/trilinos/epetra/EpetraVectorData.h"
#include "AMP/vectors/trilinos/epetra/EpetraVectorOperations.h"
#include "EpetraVector.h"


namespace AMP {
namespace LinearAlgebra {


/** \class EpetraVectorEngine
 * \brief A linear algebra engine that uses Epetra
 * \details  Use the Epetra implementation of the L1 BLAS routines.  Unlike other
 * libraries, it is very difficult to separate the data from the engine.  For this
 * reason, the EpetraVectorEngine contains the Epetra_Vector to operate on.
 */
class EpetraVectorEngine : public Vector
{


public:
    /** \brief Constructor
     * \param[in]  alias  The parameters to construct this engine
     * \param[in]  p  The buffer to use to construct the engine
     */
    explicit EpetraVectorEngine( std::shared_ptr<EpetraVectorEngineParameters> alias,
                                 std::shared_ptr<VectorData> p = nullptr );

    /** \brief Destructor
     */
    virtual ~EpetraVectorEngine() {}
#if 1
    /** \brief  Get the raw Epetra_Vector
     * \return  The Epetra_Vector currently used by this engine
     */
    inline Epetra_Vector &getEpetra_Vector() { return dynamic_cast<EpetraVectorData *>(d_VectorData)->getEpetra_Vector(); }

    /** \brief  Get the raw Epetra_Vector
     * \return  The Epetra_Vector currently used by this engine
     */
    inline const Epetra_Vector &getEpetra_Vector() const { return dynamic_cast<EpetraVectorData *>(d_VectorData)->getEpetra_Vector(); }
#else
    /** \brief  Get the raw Epetra_Vector
     * \return  The Epetra_Vector currently used by this engine
     */
    inline Epetra_Vector &getEpetra_Vector() { return d_epetraVector; }

    /** \brief  Get the raw Epetra_Vector
     * \return  The Epetra_Vector currently used by this engine
     */
    inline const Epetra_Vector &getEpetra_Vector() const { return d_epetraVector; }
#endif
public: // Functions derived from VectorData
    using VectorData::getComm;
    AMP_MPI getComm() const override { return d_Params->getComm(); }

public: // Functions derived from Vector
    using Vector::cloneVector;

    std::string type() const override { return "EpetraVectorEngine"; }
    Vector::shared_ptr cloneVector( const Variable::shared_ptr name ) const override;
    void swapVectors( Vector &other ) override;
    void aliasVector( Vector &other ) override;
    void assemble() override;

 /****************************************************************
 * VectorData operations -- will move to Vector eventually      *
 ****************************************************************/
public:
    std::string VectorDataName() const override { return d_VectorData->VectorDataName(); }
    size_t numberOfDataBlocks() const override { return d_VectorData->numberOfDataBlocks(); }
    size_t sizeOfDataBlock( size_t i = 0 ) const override { return d_VectorData->sizeOfDataBlock(i); }
    void putRawData( const double *buf ) override { d_VectorData->putRawData(buf); }
    void copyOutRawData( double *buf ) const override { d_VectorData->copyOutRawData(buf); } 
    size_t getLocalSize() const override { return d_VectorData->getLocalSize(); } 
    size_t getGlobalSize() const override { return d_VectorData->getGlobalSize(); } 
    size_t getLocalStartID() const override { return d_VectorData->getLocalStartID(); } 
    void setValuesByLocalID( int num, size_t *indices, const double *vals ) override { d_VectorData->setValuesByLocalID(num, indices, vals); }
    void setLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) override { d_VectorData->setLocalValuesByGlobalID(num, indices, vals); }
    void addValuesByLocalID( int num, size_t *indices, const double *vals ) override { d_VectorData->addValuesByLocalID(num, indices, vals); }
    void addLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) override { d_VectorData->addLocalValuesByGlobalID(num, indices, vals); }
    void getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const override { d_VectorData->getLocalValuesByGlobalID(num, indices, vals); }
    uint64_t getDataID() const override { return d_VectorData->getDataID(); }
    void *getRawDataBlockAsVoid( size_t i ) override { return d_VectorData->getRawDataBlockAsVoid(i); }
    const void *getRawDataBlockAsVoid( size_t i ) const override { return d_VectorData->getRawDataBlockAsVoid(i); }
    size_t sizeofDataBlockType( size_t i ) const override { return d_VectorData->sizeofDataBlockType(i); }
    bool isTypeId( size_t hash, size_t block ) const override { return d_VectorData->isTypeId(hash, block); }
    void swapData( VectorData &rhs ) override { d_VectorData->swapData(rhs); }
    std::shared_ptr<VectorData> cloneData() const override { return d_VectorData->cloneData(); }
    UpdateState getUpdateStatus() const override { return d_VectorData->getUpdateStatus(); }
    void setUpdateStatus( UpdateState state ) override { d_VectorData->setUpdateStatus(state); }
    void makeConsistent( ScatterType t ) override { d_VectorData->makeConsistent(t); }
    //    AMP_MPI getComm() const override{ return d_VectorData->getComm(); }
    bool hasComm() const override { return d_VectorData->hasComm(); }
    //! Get the CommunicationList for this Vector
    CommunicationList::shared_ptr getCommunicationList() const override { return d_VectorData->getCommunicationList(); }
    void setCommunicationList( CommunicationList::shared_ptr comm ) override { d_VectorData->setCommunicationList(comm); }
    void dataChanged() override { return d_VectorData->dataChanged(); }

    // missed on first round
    size_t getGlobalMaxID() const override { return d_VectorData->getGlobalMaxID(); }
    size_t getLocalMaxID() const override { return d_VectorData->getLocalMaxID(); }
    size_t getGhostSize() const override { return d_VectorData->getGhostSize(); }
    void setGhostValuesByGlobalID( int num, size_t *indices, const double *vals ) override { d_VectorData->setGhostValuesByGlobalID(num, indices, vals); }
    void setValuesByGlobalID( int num, size_t *indices, const double *vals ) override { d_VectorData->setValuesByGlobalID(num, indices, vals); }
    void addValuesByGlobalID( int num, size_t *indices, const double *vals ) override { d_VectorData->addValuesByGlobalID(num, indices, vals); }
    void getGhostAddValuesByGlobalID( int num, size_t *indices, double *vals ) const override { d_VectorData->getGhostAddValuesByGlobalID(num, indices, vals); }
    void getValuesByGlobalID( int num, size_t *indices, double *vals ) const override { d_VectorData->getValuesByGlobalID(num, indices, vals); }
    void getGhostValuesByGlobalID( int num, size_t *indices, double *vals ) const override { d_VectorData->getGhostValuesByGlobalID(num, indices, vals); }
    void getValuesByLocalID( int num, size_t *indices, double *vals ) const override { d_VectorData->getValuesByLocalID(num, indices, vals); }
/****************************************************************
 ****************************************************************/

protected:
    std::shared_ptr<EpetraVectorEngineParameters> d_Params;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
