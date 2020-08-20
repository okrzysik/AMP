#ifndef included_AMP_ArrayVector
#define included_AMP_ArrayVector

#include <string>

#include "AMP/utils/Array.h"
#include "AMP/utils/FunctionTable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/ArrayVectorData.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"
#include "AMP/vectors/operations/VectorOperationsDefault.hpp"


namespace AMP {

namespace LinearAlgebra {

/** \brief A core-local vector
 * \details This is a Vector that implements the Vector interface for a std::vector<double>.
 */
template<typename T, typename FUN = FunctionTable, typename Allocator = std::allocator<T>>
class ArrayVector : public Vector
{
private:

    ArrayVector();
    ArrayVector( const ArrayVector & );

public:
    ArrayVector( std::shared_ptr<ArrayVectorData<T, FUN, Allocator>> data );

    /** \brief    Create a ArrayVector
     * \details  This is the factory method for the ArrayVector.  It returns the shared pointer
     * to be used in the code
     * \param    localSize  The number of elements in the vector on this processor
     * \param    var The variable associated with the new vector
     */
    static Vector::shared_ptr create( const std::vector<size_t> &localSize,
                                      Variable::shared_ptr var );

    /** \brief    Create a ArrayVector
     * \details  This is the factory method for the ArrayVector.  It returns the shared pointer
     * to be used in the code
     * \param    localSize  The number of elements in the vector on this processor
     * \param    var The variable associated with the new vector
     * \param    comm The variable associated with the new vector
     */
    static Vector::shared_ptr
    create( const std::vector<size_t> &localSize, Variable::shared_ptr var, AMP_MPI comm );

    /** \brief    Create a ArrayVector
     * \details  This is the factory method for the ArrayVector.  It returns the shared pointer
     * to be used in the code that spans a comm and contains ghost values.
     * \param    var The variable associated with the new vector
     * \param    DOFs The DOFManager
     * \param    commlist The communication list
     */
    static Vector::shared_ptr create( Variable::shared_ptr var,
                                      AMP::Discretization::DOFManager::shared_ptr DOFs,
                                      AMP::LinearAlgebra::CommunicationList::shared_ptr commlist );

    /** \brief  Destructor
     */
    virtual ~ArrayVector() {}

    std::string type() const override { return "ArrayVector"; }

    using Vector::cloneVector;
    Vector::shared_ptr cloneVector( const Variable::shared_ptr name ) const override;
    void swapVectors( Vector &other ) override;
    void aliasVector( Vector &other ) override;
    /**
     * \brief This method is used to implement the assemble interface
     * of PETSc.
     * \details  This method is empty except for instantiations of NativePetscVector
     */
    void assemble() override { AMP_ERROR( "Not implemented" ); }

    //! resize the ArrayVector and reset the internal data structures
    void resize( const std::vector<size_t> &localDims );

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
};

} // namespace LinearAlgebra
} // namespace AMP

#include "ArrayVector.hpp"

#endif
