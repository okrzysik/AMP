#ifndef included_AMP_NativePetscVector
#define included_AMP_NativePetscVector

#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/petsc/PetscVector.h"
#include "AMP/vectors/petsc/NativePetscVectorData.h"

namespace AMP {
namespace LinearAlgebra {


/*! \struct Vec
    \brief PETSc vector
*/


/** \class NativePetscVector
 * \brief An AMP Vector that uses PETSc for parallel data management, linear algebra,
 * etc.
 * \details  This is an AMP wrapper to PETSc.  This is different from ManagedPetscVector
 * in that this class does not replace calls to Vec*.  Rather, it wraps these calls.
 * This class is used when PETSc is chosen as the default linear algebra engine.
 *
 * This class is not to be used directly, just through base class interfaces.
 * \see PetscVector
 * \see ManagedPetscVector
 */
class NativePetscVector : public Vector
{
public:
    /** \brief Construct a wrapper for a PETSc Vec from a set of parameters
     * \param[in] params The parameters describing the Vec
     */
    explicit NativePetscVector( VectorParameters::shared_ptr params );
    explicit NativePetscVector( std::shared_ptr<VectorData> data );

    //! Destructor
    virtual ~NativePetscVector();

    std::string type() const override { return "Native PETSc Vector"; }

    using Vector::cloneVector;
    Vector::shared_ptr cloneVector( const Variable::shared_ptr ) const override;
    void swapVectors( Vector &other ) override;
    void aliasVector( Vector & ) override;
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
/****************************************************************
 ****************************************************************/
#if 0
    std::string VectorDataName() const override { return "NativePetscVector"; }
    size_t numberOfDataBlocks() const override;
    size_t sizeOfDataBlock( size_t i ) const override;
    void putRawData( const double * ) override;
    void copyOutRawData( double *out ) const override;
    size_t getLocalSize() const override;
    size_t getGlobalSize() const override;

    void setValuesByLocalID( int, size_t *, const double * ) override;
    void setLocalValuesByGlobalID( int, size_t *, const double * ) override;
    void addValuesByLocalID( int, size_t *, const double * ) override;
    void addLocalValuesByGlobalID( int, size_t *, const double * ) override;

    void getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const override;
    void getValuesByLocalID( int numVals, size_t *ndx, double *vals ) const override;

    // Return the id of the data
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
    void swapData( VectorData &rhs ) override;

#endif



    std::shared_ptr<ParameterBase> getParameters() override;

private:
    std::shared_ptr<VectorData> d_VectorDataSP; // shared pointer to maintain scope, will go away hopefully

#if 0
    // We can always delete a NativePetscVector
    bool petscHoldsView() const override { return false; }

protected:

    void resetArray();
    void resetArray() const;

private:
    friend class NativePetscVectorOperations;
    std::shared_ptr<VectorParameters> d_pParameters;
    bool d_bDeleteMe;
    mutable double *d_pArray; // mutable so that we can cache the value
#endif
};


} // namespace LinearAlgebra
} // namespace AMP

#endif
