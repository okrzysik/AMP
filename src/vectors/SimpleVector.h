#ifndef included_AMP_SimpleVector
#define included_AMP_SimpleVector

#include "AMP/vectors/data/VectorDataCPU.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"
#include "AMP/vectors/operations/VectorOperationsDefault.hpp"
#include "Vector.h"


namespace AMP {
namespace LinearAlgebra {


/** \brief A core-local vector
 * \details This is a native AMP vector
 */
template<typename TYPE,
         typename VecOps  = VectorOperationsDefault<TYPE>,
         typename VecData = VectorDataCPU<TYPE>>
class SimpleVector : public Vector
{
protected:
    AMP_MPI d_comm;

    SimpleVector();
    explicit SimpleVector( const SimpleVector & );

public:
    /** \brief    Create a SimpleVector
     * \details  This is the factory method for the SimpleVector.  It returns the shared pointer
     * to be used in the code
     * \param    localSize  The number of elements in the vector on this processor
     * \param    var The variable associated with the new vector
     */
    static Vector::shared_ptr create( size_t localSize, const std::string &var );

    /** \brief    Create a SimpleVector
     * \details  This is the factory method for the SimpleVector.  It returns the shared pointer
     * to be used in the code
     * \param    localSize  The number of elements in the vector on this processor
     * \param    var The variable associated with the new vector
     */
    static Vector::shared_ptr create( size_t localSize, Variable::shared_ptr var );

    /** \brief    Create a SimpleVector
     * \details  This is the factory method for the SimpleVector.  It returns the shared pointer
     * to be used in the code
     * \param    localSize  The number of elements in the vector on this processor
     * \param    var The variable associated with the new vector
     * \param    comm The variable associated with the new vector
     */
    static Vector::shared_ptr create( size_t localSize, Variable::shared_ptr var, AMP_MPI comm );

    /** \brief    Create a SimpleVector
     * \details  This is the factory method for the SimpleVector.  It returns the shared pointer
     * to be used in the code that spans a comm and contains ghost values.
     * \param    var The variable associated with the new vector
     * \param    DOFs The DOFManager
     * \param    commlist The communication list
     */
    static Vector::shared_ptr create( Variable::shared_ptr var,
                                      AMP::Discretization::DOFManager::shared_ptr DOFs,
                                      AMP::LinearAlgebra::CommunicationList::shared_ptr commlist );

    //! Destructor
    virtual ~SimpleVector() override {}


    //! Resize this vector
    virtual void resize( size_t i );


public: // Functions derived from Vector
    using Vector::cloneVector;

    std::string type() const override { return "Simple Vector"; }
    Vector::shared_ptr cloneVector( const Variable::shared_ptr name ) const override;
    void swapVectors( Vector &other ) override;
    void aliasVector( Vector &other ) override;
    void assemble() override;

public:
/****************************************************************
 * VectorData operations -- will move to Vector eventually      *
 ****************************************************************/
    std::string VectorDataName() const override;
    size_t numberOfDataBlocks() const override;
    size_t sizeOfDataBlock( size_t i = 0 ) const override;
    void putRawData( const double *buf ) override;
    void copyOutRawData( double *buf ) const override;
    size_t getLocalSize() const override;
    size_t getGlobalSize() const override;
    size_t getLocalStartID() const override;
    void setValuesByLocalID( int num, size_t *indices, const double *vals ) override;
    void setLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) override;
    void addValuesByLocalID( int num, size_t *indices, const double *vals ) override;
    void addLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) override;
    void getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const override;
    uint64_t getDataID() const override;
    void *getRawDataBlockAsVoid( size_t i ) override;
    const void *getRawDataBlockAsVoid( size_t i ) const override;
    size_t sizeofDataBlockType( size_t i ) const override;
    bool isTypeId( size_t hash, size_t block ) const override;
    void swapData( VectorData &rhs ) override;
    std::shared_ptr<VectorData> cloneData() const override;
// temporary functions
    void setComm( AMP_MPI comm) { d_comm = comm; }
    void allocateVectorData( size_t localSize, size_t numLocal, size_t numGlobal );
void setDOFManager( std::shared_ptr<AMP::Discretization::DOFManager> dofManager);
    void setCommunicationList( CommunicationList::shared_ptr comm ) override;
    AMP_MPI Vector::getComm() const override{ return d_VectorData->getComm(); }
};


} // namespace LinearAlgebra
} // namespace AMP


#include "AMP/vectors/SimpleVector.hpp"

#endif
