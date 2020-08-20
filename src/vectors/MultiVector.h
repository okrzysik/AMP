#ifndef included_AMP_MultiVector
#define included_AMP_MultiVector


#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/MultiVectorData.h"
#include "AMP/vectors/operations/MultiVectorOperations.h"


namespace AMP {
namespace LinearAlgebra {


/** \brief A collection of AMP Vectors that appear as one vector
 * \details
 *    Given a set of vectors, they can be collected into a singlevector.  This class
 *    accomplishes this task.
 */
class MultiVector final : public Vector
{
public:
    /** \brief Return the first vector in the MultiVector
     * \return The first vector
     */
    inline auto beginVector() { return d_vVectors.begin(); }

    /** \brief Return one past the last vector in the MultiVector
     * \return One past the last vector
     */
    inline auto endVector() { return d_vVectors.end(); }

    /** \brief Determine if a Vector is a constituent
     * \param[in]  p  The vector to look for
     * \return True if the pointer p is in the multivector
     */
    bool containsPointer( const Vector::shared_ptr p ) const;

    /** \brief Create a new multivector in parallel
     * \param[in] name  Variable describing the new vector
     * \param[in] comm  Communicator to build the MultiVector on
     * \param[in] vecs  Optional list of vectors in the MultiVector
     */
    static std::shared_ptr<MultiVector>
    create( Variable::shared_ptr name,
            const AMP_MPI &comm,
            const std::vector<Vector::shared_ptr> &vecs = std::vector<Vector::shared_ptr>() );

    /** \brief Create a new multivector in parallel
     * \param[in] name  Name of the new vector
     * \param[in] comm  Communicator to build the MultiVector on
     * \param[in] vecs  Optional list of vectors in the MultiVector
     */
    static std::shared_ptr<MultiVector>
    create( const std::string &name,
            const AMP_MPI &comm,
            const std::vector<Vector::shared_ptr> &vecs = std::vector<Vector::shared_ptr>() );

    /** \brief Create a new multivector in parallel
     * \param[in] name  Variable describing the new vector
     * \param[in] comm  Communicator to build the MultiVector on
     * \param[in] vecs  Optional list of vectors in the MultiVector
     */
    static std::shared_ptr<const MultiVector>
    const_create( Variable::shared_ptr name,
                  const AMP_MPI &comm,
                  const std::vector<Vector::const_shared_ptr> &vecs =
                      std::vector<Vector::const_shared_ptr>() );

    /** \brief Create a new multivector in parallel
     * \param[in] name  Name of the new vector
     * \param[in] comm  Communicator to build the MultiVector on
     * \param[in] vecs  Optional list of vectors in the MultiVector
     */
    static std::shared_ptr<const MultiVector>
    const_create( const std::string &name,
                  const AMP_MPI &comm,
                  const std::vector<Vector::const_shared_ptr> &vecs =
                      std::vector<Vector::const_shared_ptr>() );

    /** \brief Create a multivector view of a vector
     * \param[in] vec  The vector to view
     * \param[in] comm  Communicator to create the MultiVector on
     * \details  If vec is a MultiVector, it is returned.  Otherwise, a MultiVector is created
     * and vec is added to it.  If vec is not a parallel vector(such as a SimpleVector), comm
     * must be specified.
     */
    static std::shared_ptr<MultiVector> view( Vector::shared_ptr vec,
                                              const AMP_MPI &comm = AMP_MPI( AMP_COMM_NULL ) );

    /** \brief Create a multivector view of a vector
     * \param[in] vec  The vector to view
     * \param[in] comm  Communicator to create the MultiVector on
     * \details  If vec is a MultiVector, it is returned.  Otherwise, a MultiVector is created
     * and vec is added to it.  If vec is not a parallel vector(such as a SimpleVector), comm
     * must be specified.
     */
    static std::shared_ptr<const MultiVector>
    constView( Vector::const_shared_ptr vec, const AMP_MPI &comm = AMP_MPI( AMP_COMM_NULL ) );

    /** \brief Encapsulate a vector in a MultiVector
     * \param[in] vec  The vector to view
     * \param[in] comm  Communicator to create the MultiVector on
     * \details  If vec is a MultiVector, it is returned.  Otherwise, a MultiVector is created
     * and vec is added to it.  If vec is not a parallel vector(such as a SimpleVector), comm
     * must be specified.
     */
    static std::shared_ptr<MultiVector>
    encapsulate( Vector::shared_ptr vec, const AMP_MPI &comm = AMP_MPI( AMP_COMM_NULL ) );

    /** \brief Replace a vector in a MultiVector
     * \details  This function will replace a given vector in the multivector with a different
     * vector.  The original and new vectors must share the same DOFManager.
     * This is a purly local operation and does not require communication, but should be called
     * by all processors that own the vectors.
     * \param[in] oldVec  The original vector to replace
     * \param[in] newVec  The new vector to use
     */
    virtual void replaceSubVector( Vector::shared_ptr oldVec, Vector::shared_ptr newVec );

    /** \brief  Add a vector to a MultiVector.
     *   Note: this is a collective operation, vec may be NULL on some processors
     * \param[in]  vec  The Vector to add to the MultiVector
     */
    virtual void addVector( Vector::shared_ptr vec );

    /** \brief  Add vector(s) to a MultiVector.
     *   Note: This is a collective operation, vec may be different sizes on different processors
     * \param[in]  vec  The Vector to add to the MultiVector
     */
    virtual void addVector( std::vector<Vector::shared_ptr> vec );

    /** \brief  Remove a vector from a MultiVector
     * \param[in]  vec  The Vector to remove from the MultiVector
     */
    virtual void eraseVector( Vector::shared_ptr vec );

    /** \brief  Return the i-th Vector in this MultiVector
     * \param[in]  i  Which vector to return
     * \return  The i-th Vector
     */
    virtual Vector::shared_ptr getVector( size_t i );

    /** \brief  Return the i-th Vector in this MultiVector
     * \param[in]  i  Which vector to return
     * \return  The i-th Vector
     */
    virtual Vector::const_shared_ptr getVector( size_t i ) const;

    /** \brief  Obtain the number of Vector objects in this MultiVector
     * \return  The number of Vector objects in this MultiVector
     */
    size_t getNumberOfSubvectors() const;

    //!  Destructor
    virtual ~MultiVector();
    std::string type() const override;

    AMP_MPI getComm() const override;

    Vector::shared_ptr subsetVectorForVariable( Variable::const_shared_ptr name ) override;
    Vector::const_shared_ptr
    constSubsetVectorForVariable( Variable::const_shared_ptr name ) const override;
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
    //    AMP_MPI getComm() const override { return d_VectorData->getComm(); }
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
    bool containsGlobalElement( size_t GID ) override { return d_VectorData->containsGlobalElement(GID); }

    void dumpOwnedData( std::ostream &out, size_t GIDoffset = 0, size_t LIDoffset = 0 ) const override
    {
      d_VectorData->dumpOwnedData(out, GIDoffset, LIDoffset);
    }
    void dumpGhostedData( std::ostream &out, size_t offset = 0 ) const override{ d_VectorData->dumpGhostedData(out, offset); }
/****************************************************************
 ****************************************************************/

protected:
    Vector::shared_ptr selectInto( const VectorSelector & ) override;
    Vector::const_shared_ptr selectInto( const VectorSelector &criterion ) const override;

    //! The communicator this multivector is on
    AMP_MPI d_Comm;
    //! Indicates if the multivector created the communicator
    bool d_CommCreated;

    //! The list of AMP Vectors that comprise this Vector
    std::vector<Vector::shared_ptr> d_vVectors;

    /** \brief A convenience method for extracting vectors from a base class
     * \param[in]  vec  The vector to extract a vector from
     * \param[in]  which  Which vector to get
     * \return     The extracted vector
     */
    const Vector::shared_ptr &getVector( const VectorData &vec, size_t which ) const;

    /** \brief A convenience method for extracting vectors from a base class
     * \param[in]  vec  The vector to extract a vector from
     * \param[in]  which  Which vector to get
     * \return     The extracted vector
     */
    Vector::shared_ptr &getVector( VectorData &vec, size_t which ) const;

    /** Constructor:  create a MultiVector with a particular variable
     * \param[in]  name  The vector to create the MultiVector from
     */
    explicit MultiVector( const std::string &name );


private:
    // Helper function to add a vector without updating the DOF manager
    void addVectorHelper( Vector::shared_ptr vec );

    // Helper function to reset the vector data
    inline void resetVectorData();

    // Helper function to reset the vector operations
    inline void resetVectorOperations();

public: // Pull Vector into the current scope
    using Vector::cloneVector;
    using Vector::constSubsetVectorForVariable;
    using Vector::subsetVectorForVariable;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
