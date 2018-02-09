#ifndef included_AMP_MultiVector
#define included_AMP_MultiVector


#include "AMP/vectors/DataChangePassThrough.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorEngine.h"
#include "AMP/vectors/data/MultiVectorData.h"
#include "AMP/vectors/operations/MultiVectorOperations.h"


namespace AMP {
namespace LinearAlgebra {


/** \brief A collection of AMP Vectors that appear as one vector
 * \details
 *    Given a set of vectors, they can be collected into a singlevector.  This class
 *    accomplishes this task.
 */
class MultiVector : public Vector,
                    public MultiVectorData,
                    public MultiVectorOperations,
                    public VectorEngine,
                    public DataChangePassThrough
{
public:
    //!  Iterator typedef
    typedef std::vector<Vector::shared_ptr>::iterator vector_iterator;

    //!  Iterator typedef
    typedef std::vector<Vector::shared_ptr>::const_iterator vector_const_iterator;

    /** \brief Return the first vector in the MultiVector
     * \return The first vector
     */
    vector_iterator beginVector();

    /** \brief Return one past the last vector in the MultiVector
     * \return One past the last vector
     */
    vector_iterator endVector();

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
    static AMP::shared_ptr<MultiVector>
    create( Variable::shared_ptr name,
            AMP_MPI comm,
            const std::vector<Vector::shared_ptr> &vecs = std::vector<Vector::shared_ptr>() );

    /** \brief Create a new multivector in parallel
     * \param[in] name  Name of the new vector
     * \param[in] comm  Communicator to build the MultiVector on
     * \param[in] vecs  Optional list of vectors in the MultiVector
     */
    static AMP::shared_ptr<MultiVector>
    create( const std::string &name,
            AMP_MPI comm,
            const std::vector<Vector::shared_ptr> &vecs = std::vector<Vector::shared_ptr>() );

    /** \brief Create a new multivector in parallel
     * \param[in] name  Variable describing the new vector
     * \param[in] comm  Communicator to build the MultiVector on
     * \param[in] vecs  Optional list of vectors in the MultiVector
     */
    static AMP::shared_ptr<const MultiVector>
    const_create( Variable::shared_ptr name,
                  AMP_MPI comm,
                  const std::vector<Vector::const_shared_ptr> &vecs =
                      std::vector<Vector::const_shared_ptr>() );

    /** \brief Create a new multivector in parallel
     * \param[in] name  Name of the new vector
     * \param[in] comm  Communicator to build the MultiVector on
     * \param[in] vecs  Optional list of vectors in the MultiVector
     */
    static AMP::shared_ptr<const MultiVector>
    const_create( const std::string &name,
                  AMP_MPI comm,
                  const std::vector<Vector::const_shared_ptr> &vecs =
                      std::vector<Vector::const_shared_ptr>() );

    /** \brief Create a multivector view of a vector
     * \param[in] vec  The vector to view
     * \param[in] comm  Communicator to create the MultiVector on
     * \details  If vec is a MultiVector, it is returned.  Otherwise, a MultiVector is created
     * and vec is added to it.  If vec is not a parallel vector(such as a SimpleVector), comm
     * must be specified.
     */
    static AMP::shared_ptr<MultiVector> view( Vector::shared_ptr vec,
                                              AMP_MPI comm = AMP_MPI( AMP_COMM_NULL ) );

    /** \brief Create a multivector view of a vector
     * \param[in] vec  The vector to view
     * \param[in] comm  Communicator to create the MultiVector on
     * \details  If vec is a MultiVector, it is returned.  Otherwise, a MultiVector is created
     * and vec is added to it.  If vec is not a parallel vector(such as a SimpleVector), comm
     * must be specified.
     */
    static AMP::shared_ptr<const MultiVector> constView( Vector::const_shared_ptr vec,
                                                         AMP_MPI comm = AMP_MPI( AMP_COMM_NULL ) );

    /** \brief Encapsulate a vector in a MultiVector
     * \param[in] vec  The vector to view
     * \param[in] comm  Communicator to create the MultiVector on
     * \details  If vec is a MultiVector, it is returned.  Otherwise, a MultiVector is created
     * and vec is added to it.  If vec is not a parallel vector(such as a SimpleVector), comm
     * must be specified.
     */
    static AMP::shared_ptr<MultiVector> encapsulate( Vector::shared_ptr vec,
                                                     AMP_MPI comm = AMP_MPI( AMP_COMM_NULL ) );

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

    virtual void
    dumpOwnedData( std::ostream &out, size_t GIDoffset = 0, size_t LIDoffset = 0 ) const override;
    virtual void dumpGhostedData( std::ostream &out, size_t offset = 0 ) const override;
    virtual std::string type() const override;

    virtual AMP_MPI getComm() const override;

    virtual Vector::shared_ptr subsetVectorForVariable( Variable::const_shared_ptr name ) override;
    virtual Vector::const_shared_ptr
    constSubsetVectorForVariable( Variable::const_shared_ptr name ) const override;
    virtual Vector::shared_ptr cloneVector( const Variable::shared_ptr name ) const override;
    virtual void swapVectors( Vector &other ) override;
    virtual void aliasVector( Vector &other ) override;
    virtual void assemble() override;

    // Vector engine functions
    virtual AMP::shared_ptr<VectorData> getNewBuffer() override;
    virtual bool sameEngine( VectorEngine &rhs ) const override;
    virtual AMP::shared_ptr<VectorEngine> cloneEngine( AMP::shared_ptr<VectorData> ) const override;
    virtual void swapEngines( AMP::shared_ptr<VectorEngine> p ) override;


protected:
    virtual Vector::shared_ptr selectInto( const VectorSelector & ) override;
    virtual Vector::const_shared_ptr selectInto( const VectorSelector &criterion ) const override;

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
    const Vector::shared_ptr &getVector( const VectorOperations &vec, size_t which ) const;

    /** \brief A convenience method for extracting vectors from a base class
     * \param[in]  vec  The vector to extract a vector from
     * \param[in]  which  Which vector to get
     * \return     The extracted vector
     */
    Vector::shared_ptr &getVector( VectorOperations &vec, size_t which ) const;

    /** Constructor:  create a MultiVector with a particular variable
     * \param[in]  name  The vector to create the MultiVector from
     */
    explicit MultiVector( const std::string &name );

    virtual void dataChanged() override;

private:
    // Helper function to add a vector without updating the DOF manager
    void addVectorHelper( Vector::shared_ptr vec );

    // Helper function to update the vector data
    inline void updateVectorData();

    // Helper function to update the vector operations
    inline void updateVectorOperations();


public: // Pull Vector into the current scope
    using Vector::cloneVector;
    using Vector::constSubsetVectorForVariable;
    using Vector::subsetVectorForVariable;


public: // Pull VectorOperations into the current scope
    using MultiVectorOperations::abs;
    using MultiVectorOperations::add;
    using MultiVectorOperations::addScalar;
    using MultiVectorOperations::axpby;
    using MultiVectorOperations::axpy;
    using MultiVectorOperations::divide;
    using MultiVectorOperations::dot;
    using MultiVectorOperations::equals;
    using MultiVectorOperations::linearSum;
    using MultiVectorOperations::minQuotient;
    using MultiVectorOperations::multiply;
    using MultiVectorOperations::reciprocal;
    using MultiVectorOperations::scale;
    using MultiVectorOperations::setRandomValues;
    using MultiVectorOperations::subtract;
    using MultiVectorOperations::wrmsNorm;
    using MultiVectorOperations::wrmsNormMask;
    using MultiVectorOperations::zero;
};


} // namespace LinearAlgebra
} // namespace AMP

#include "MultiVector.inline.h"


#endif
