#ifndef included_AMP_NativePetscVector
#define included_AMP_NativePetscVector

#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"
#include "AMP/vectors/petsc/PetscVector.h"

namespace AMP {
namespace LinearAlgebra {


/*! \struct Vec
    \brief PETSc vector
*/


/** \class NativePetscVectorParameters
 * \brief Parameters to set when creating a NativePetscVector
 */
class NativePetscVectorParameters : public VectorParameters
{
public:
    //!  The vector to wrap
    Vec d_InVec;

    //!  The communicator associated with the Vec
    AMP_MPI d_Comm;

    //!  If true, ~NativePetscVector() will call VecDestroy()
    bool d_Deleteable;

    //! The number of local entities in the vector
    size_t d_localsize;

    /** \brief Constructor
     * \details  This will create NativePetscVectorParameters that can be used
     *   to construct a NativePetscVector around the given ::Vec.  Note that the
     *   existing vector must be destroyed once through a call to VecDestroy.
     *   This can be done be the user or the NativePetscVector based on the deleteable
     *   flag.
     * \param[in] v            The vector to wrap
     * \param[in] deleteable   Do we want ~NativePetscVector() to call VecDestroy() on v
     */
    NativePetscVectorParameters( Vec v, bool deleteable );
};


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
class NativePetscVector : public Vector,
                          public PetscVector,
                          public VectorOperationsDefault<double>
{
public:

    /** \brief Construct a wrapper for a PETSc Vec from a set of parameters
     * \param[in] params The parameters describing the Vec
     */
    explicit NativePetscVector( VectorParameters::shared_ptr params );

    //! Destructor
    virtual ~NativePetscVector();

    std::string type() const override { return "Native PETSc Vector"; }
    
    using Vector::cloneVector;
    Vector::shared_ptr cloneVector( const Variable::shared_ptr ) const override;
    void copy( const VectorOperations &vec ) override;

    void swapVectors( Vector &other ) override;
    void aliasVector( Vector & ) override;

    std::string VectorDataName() const override { return "NativePetscVector"; }
    size_t numberOfDataBlocks() const override;
    size_t sizeOfDataBlock( size_t i ) const override;


    void setToScalar( double alpha ) override;
    void scale( double alpha, const VectorOperations &x ) override;
    void scale( double alpha ) override;
    void add( const VectorOperations &x, const VectorOperations &y ) override;
    void subtract( const VectorOperations &x, const VectorOperations &y ) override;
    void multiply( const VectorOperations &x, const VectorOperations &y ) override;
    void divide( const VectorOperations &x, const VectorOperations &y ) override;
    void reciprocal( const VectorOperations &x ) override;
    void linearSum( double alpha,
                    const VectorOperations &x,
                    double beta,
                    const VectorOperations &y ) override;
    void axpy( double alpha, const VectorOperations &x, const VectorOperations &y ) override;
    void axpby( double alpha, double beta, const VectorOperations &x ) override;
    void abs( const VectorOperations &x ) override;
    double min( void ) const override;
    double max( void ) const override;
    void setRandomValues( void ) override;

    double L1Norm( void ) const override;
    double L2Norm( void ) const override;
    double maxNorm( void ) const override;
    using Vector::dot;
    double dot( const VectorOperations &x ) const override;

    double localL1Norm( void ) const override;
    double localL2Norm( void ) const override;
    double localMaxNorm( void ) const override;

    void setValuesByLocalID( int, size_t *, const double * ) override;
    void setLocalValuesByGlobalID( int, size_t *, const double * ) override;
    void addValuesByLocalID( int, size_t *, const double * ) override;
    void addLocalValuesByGlobalID( int, size_t *, const double * ) override;

    void getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const override;
    void getValuesByLocalID( int numVals, size_t *ndx, double *vals ) const override;

    void assemble() override;

    size_t getLocalSize() const override;
    size_t getGlobalSize() const override;

    void putRawData( const double * ) override;
    void swapData( VectorData &rhs ) override;

    void copyOutRawData( double *out ) const override;

    std::shared_ptr<ParameterBase> getParameters() override;

    // We can always delete a NativePetscVector
    bool petscHoldsView() const override { return false; }

    // Return the id of the data
    uint64_t getDataID() const override
    {
        return reinterpret_cast<uint64_t>( getRawDataBlockAsVoid( 0 ) );
    }

    bool isTypeId( size_t hash, size_t ) const override
    {
        return hash == typeid( double ).hash_code();
    }
    size_t sizeofDataBlockType( size_t ) const override { return sizeof( double ); }

protected:
    void *getRawDataBlockAsVoid( size_t i ) override;
    const void *getRawDataBlockAsVoid( size_t i ) const override;

    void resetArray();
    void resetArray() const;

    // Function to perform  this = alpha x + beta y + gamma this
    virtual void axpbypcz( double alpha,
                           const VectorOperations &x,
                           double beta,
                           const VectorOperations &y,
                           double gamma );

private:
    std::shared_ptr<VectorParameters> d_pParameters;
    bool d_bDeleteMe;
    mutable double *d_pArray; // mutable so that we can cache the value

public: // Pull VectorOperations into the current scope
    using VectorOperationsDefault::minQuotient;
    using VectorOperationsDefault::wrmsNorm;
    using VectorOperationsDefault::wrmsNormMask;
};


} // namespace LinearAlgebra
} // namespace AMP

#endif
