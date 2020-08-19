#ifndef included_AMP_ManagedVector
#define included_AMP_ManagedVector

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/ManagedVectorData.h"

#include <stdexcept>
#include <vector>


namespace AMP {
namespace LinearAlgebra {

/**
   \brief Class used to control data and kernels of various vector libraries
   \details  A ManagedVector will take an engine and create a buffer, if
   necessary.

   A ManagedVector has two pointers: data and engine.  If the data pointer
   is null, then the engine is assumed to have the data.
*/
class ManagedVector : public Vector
{

public:
    /** \brief Construct a ManagedVector from a set of parameters
     * \param[in] params  The description of the ManagedVector
     */
    explicit ManagedVector( VectorParameters::shared_ptr params );

    /** \brief Construct a view of an AMP vector
     * \param[in] alias  Vector to view
     */
    explicit ManagedVector( const Vector::shared_ptr alias );

    //! Destructor
    virtual ~ManagedVector();

    /** \brief  If a vector has multiple views to multiple external packages
     * associated with it, this will return the barest version of the vector
     * \return A vector with the fewest views associated with it.
     * \details  A ManagedVector must have an engine and it may have data.
     * If it has an engine with no data, then the engine has must have data.
     * If the engine can be cast to a ManagedVector, it is and getRootVector
     * is called recursively.
     */
    Vector::shared_ptr getRootVector();

    /** \brief  Return the engine associated with this ManagedVector
     * \return The engine
     */
    std::shared_ptr<Vector> getVectorEngine();
    std::shared_ptr<const Vector> getVectorEngine() const;

    virtual bool isAnAliasOf( Vector &rhs );
    virtual bool isAnAliasOf( Vector::shared_ptr rhs );

    virtual std::shared_ptr<ManagedVectorParameters> getManagedVectorParameters();

protected:

    //! The parameters used to create this vector
    std::shared_ptr<ManagedVectorParameters> d_pParameters;

    //! Function that returns a pointer to a managed vector
    virtual ManagedVector *getNewRawPtr() const = 0;

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
    AMP_MPI getComm() const override { return d_VectorData->getComm(); }
    //! Get the CommunicationList for this Vector
    CommunicationList::shared_ptr getCommunicationList() const override { return d_VectorData->getCommunicationList(); }
    void setCommunicationList( CommunicationList::shared_ptr comm ) override { d_VectorData->setCommunicationList(comm); }
    void dataChanged() override { return d_VectorData->dataChanged(); }
/****************************************************************
 ****************************************************************/

public: // Derived from Vector
    using Vector::cloneVector;
    std::string type() const override;
    std::shared_ptr<Vector> cloneVector( const Variable::shared_ptr name ) const override;
    std::shared_ptr<ParameterBase> getParameters() override;
    Vector::shared_ptr subsetVectorForVariable( Variable::const_shared_ptr name ) override;
    Vector::const_shared_ptr
    constSubsetVectorForVariable( Variable::const_shared_ptr name ) const override;
    void swapVectors( Vector &other ) override;
    void aliasVector( Vector &other ) override;

protected: // Derived from Vector
    Vector::shared_ptr selectInto( const VectorSelector & ) override;
    Vector::const_shared_ptr selectInto( const VectorSelector & ) const override;

private:
    ManagedVector();
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
