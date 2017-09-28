#ifndef included_AMP_NativePetscVector
#define included_AMP_NativePetscVector

#include "utils/AMP_MPI.h"
#include "vectors/NativeVector.h"
#include "vectors/VectorEngine.h"
#include "vectors/operations/VectorOperationsDefault.h"
#include "vectors/petsc/PetscVector.h"


namespace AMP {
namespace LinearAlgebra {


/*! \struct Vec
    \brief PETSc vector
*/


/** \class NativePetscVectorParameters
 * \brief Parameters to set when creating a NativePetscVector
 */
class NativePetscVectorParameters : public NativeVectorParameters
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
class NativePetscVector : public NativeVector,
                          public PetscVector,
                          public VectorEngine,
                          public VectorOperationsDefault<double>
{
public:
    //! Conveninece typedef
    typedef NativeVector::parameters_ptr parameters_ptr;

    //! Conveninece typedef
    typedef NativeVectorParameters parameters;


    /** \brief Construct a wrapper for a PETSc Vec from a set of parameters
     * \param[in] params The parameters describing the Vec
     */
    explicit NativePetscVector( VectorParameters::shared_ptr params );

    //! Destructor
    virtual ~NativePetscVector();

    virtual std::string type() const override { return "Native PETSc Vector"; }

    virtual Vector::shared_ptr getManagedVectorCopy( AMP_MPI comm ) override;

    virtual Vector::shared_ptr getManagedVectorDuplicate( AMP_MPI comm ) override;

    using Vector::cloneVector;
    virtual Vector::shared_ptr cloneVector( const Variable::shared_ptr ) const override;
    virtual void copy( const VectorOperations &vec ) override;

    virtual void swapVectors( Vector &other ) override;
    virtual void aliasVector( Vector & ) override;

    virtual std::string VectorDataName() const override { return "NativePetscVector"; }
    virtual size_t numberOfDataBlocks() const override;
    virtual size_t sizeOfDataBlock( size_t i ) const override;


    virtual void setToScalar( double alpha ) override;
    virtual void scale( double alpha, const VectorOperations &x ) override;
    virtual void scale( double alpha ) override;
    virtual void add( const VectorOperations &x, const VectorOperations &y ) override;
    virtual void subtract( const VectorOperations &x, const VectorOperations &y ) override;
    virtual void multiply( const VectorOperations &x, const VectorOperations &y ) override;
    virtual void divide( const VectorOperations &x, const VectorOperations &y ) override;
    virtual void reciprocal( const VectorOperations &x ) override;
    virtual void linearSum( double alpha,
                            const VectorOperations &x,
                            double beta,
                            const VectorOperations &y ) override;
    virtual void
    axpy( double alpha, const VectorOperations &x, const VectorOperations &y ) override;
    virtual void axpby( double alpha, double beta, const VectorOperations &x ) override;
    virtual void abs( const VectorOperations &x ) override;
    virtual double min( void ) const override;
    virtual double max( void ) const override;
    virtual void setRandomValues( void ) override;

    virtual double L1Norm( void ) const override;
    virtual double L2Norm( void ) const override;
    virtual double maxNorm( void ) const override;
    using Vector::dot;
    virtual double dot( const VectorOperations &x ) const override;

    virtual double localL1Norm( void ) const override;
    virtual double localL2Norm( void ) const override;
    virtual double localMaxNorm( void ) const override;

    virtual void setValuesByLocalID( int, size_t *, const double * ) override;
    virtual void setLocalValuesByGlobalID( int, size_t *, const double * ) override;
    virtual void addValuesByLocalID( int, size_t *, const double * ) override;
    virtual void addLocalValuesByGlobalID( int, size_t *, const double * ) override;

    virtual void getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const override;
    virtual void getValuesByLocalID( int numVals, size_t *ndx, double *vals ) const override;

    virtual void assemble() override;

    virtual size_t getLocalSize() const override;
    virtual size_t getGlobalSize() const override;

    virtual void putRawData( const double * ) override;

    virtual AMP::shared_ptr<VectorData> getNewBuffer() override;
    virtual bool sameEngine( VectorEngine & ) const override;
    virtual AMP::shared_ptr<VectorEngine>
    cloneEngine( AMP::shared_ptr<VectorData> p ) const override;

    virtual void swapEngines( AMP::shared_ptr<VectorEngine> ) override;
    virtual void swapData( VectorData &rhs ) override;

    virtual AMP_MPI getComm() const override;

    virtual void copyOutRawData( double *out ) const override;

    virtual AMP::shared_ptr<ParameterBase> getParameters() override;

    // We can always delete a NativePetscVector
    virtual bool petscHoldsView() const override { return false; }

    // Return the id of the data
    virtual uint64_t getDataID() const override
    {
        return reinterpret_cast<uint64_t>( getRawDataBlockAsVoid( 0 ) );
    }

    virtual bool isTypeId( size_t hash, size_t ) const override
    {
        return hash == typeid( double ).hash_code();
    }
    virtual size_t sizeofDataBlockType( size_t ) const override { return sizeof( double ); }

protected:
    virtual void *getRawDataBlockAsVoid( size_t i ) override;
    virtual const void *getRawDataBlockAsVoid( size_t i ) const override;

    void resetArray();
    void resetArray() const;

    // Function to perform  this = alpha x + beta y + gamma this
    virtual void axpbypcz( double alpha,
                           const VectorOperations &x,
                           double beta,
                           const VectorOperations &y,
                           double gamma );

private:
    parameters_ptr d_pParameters;
    bool d_bDeleteMe;
    mutable double *d_pArray; // mutable so that we can cache the value

public: // Pull VectorOperations into the current scope
    using VectorOperationsDefault::abs;
    using VectorOperationsDefault::add;
    using VectorOperationsDefault::axpby;
    using VectorOperationsDefault::axpy;
    using VectorOperationsDefault::divide;
    using VectorOperationsDefault::dot;
    using VectorOperationsDefault::linearSum;
    using VectorOperationsDefault::minQuotient;
    using VectorOperationsDefault::multiply;
    using VectorOperationsDefault::reciprocal;
    using VectorOperationsDefault::scale;
    using VectorOperationsDefault::setRandomValues;
    using VectorOperationsDefault::subtract;
    using VectorOperationsDefault::wrmsNorm;
    using VectorOperationsDefault::wrmsNormMask;
};


} // namespace LinearAlgebra
} // namespace AMP

#include "NativePetscVector.inline.h"

#endif
