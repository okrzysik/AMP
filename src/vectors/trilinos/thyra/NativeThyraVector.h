#ifndef included_AMP_NativeThyraVector
#define included_AMP_NativeThyraVector

#include "AMP/vectors/NativeVector.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"
#include "AMP/vectors/trilinos/thyra/ThyraVector.h"

namespace AMP {
namespace LinearAlgebra {


/** \class NativeThyraVectorParameters
 * \brief Parameters to set when creating a NativeThyraVector
 */
class NativeThyraVectorParameters : public NativeVectorParameters
{
public:
    //! The vector to wrap
    Teuchos::RCP<Thyra::VectorBase<double>> d_InVec;

    //! The local size of the vector
    size_t d_local;

    //! The comm of the vector
    AMP_MPI d_comm;

    //! The variable to use with the vector
    Variable::shared_ptr d_var;
};


/** \class NativeThyraVector
 * \brief An AMP Vector that uses Thyra for parallel data management, linear algebra,
 * etc.
 * \details  This is an AMP wrapper to Thyra.  This is different from ManagedThyraVector
 * in that this class does not replace calls to Vec*.  Rather, it wraps these calls.
 * This class is used when Thyra is chosen as the default linear algebra engine.
 *
 * This class is not to be used directly, just through base class interfaces.
 * \see ThyraVector
 * \see ManagedThyraVector
 */
class NativeThyraVector : public NativeVector,
                          public ThyraVector,
                          public VectorOperationsDefault<double>
{
public:
    //! Conveninece typedef
    typedef NativeVector::parameters_ptr parameters_ptr;

    //! Conveninece typedef
    typedef NativeVectorParameters parameters;


    /** \brief Construct a wrapper for a Thyra Vec from a set of parameters
     * \param[in] params The parameters describing the Vec
     */
    explicit NativeThyraVector( VectorParameters::shared_ptr params );

    //! Destructor
    virtual ~NativeThyraVector();


    //! Overloaded functions
    std::string type() const override { return "Native Thyra Vector"; }
    std::string VectorDataName() const override { return "NativeThyraVector"; }
    Vector::shared_ptr getManagedVectorCopy( AMP_MPI comm ) override;
    Vector::shared_ptr getManagedVectorDuplicate( AMP_MPI comm ) override;
    Vector::shared_ptr cloneVector( const Variable::shared_ptr ) const override;
    void copy( const VectorOperations &vec ) override;
    void swapVectors( Vector &other ) override;
    void aliasVector( Vector & ) override;
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
    void
    axpy( double alpha, const VectorOperations &x, const VectorOperations &y ) override;
    void axpby( double alpha, double beta, const VectorOperations &x ) override;
    void abs( const VectorOperations &x ) override;
    double min( void ) const override;
    double max( void ) const override;
    void setRandomValues( void ) override;
    double L1Norm( void ) const override;
    double L2Norm( void ) const override;
    double maxNorm( void ) const override;
    double dot( const VectorOperations &x ) const override;
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
    void copyOutRawData( double *out ) const override;
    uint64_t getDataID() const override
    {
        return reinterpret_cast<uint64_t>( getRawDataBlockAsVoid( 0 ) );
    }
    bool isTypeId( size_t hash, size_t ) const override
    {
        return hash == typeid( double ).hash_code();
    }
    void swapData( VectorData & ) override { AMP_ERROR( "Not finished" ); }


protected:
    //! Empty constructor.
    NativeThyraVector();

    void *getRawDataBlockAsVoid( size_t i ) override;
    const void *getRawDataBlockAsVoid( size_t i ) const override;
    size_t sizeofDataBlockType( size_t ) const override { return sizeof( double ); }

private:
    size_t d_local;

    static Teuchos::RCP<const Thyra::VectorBase<double>> getThyraVec( const VectorOperations &v );
    static Teuchos::RCP<const Thyra::VectorBase<double>>
    getThyraVec( const Vector::const_shared_ptr &v );


public: // Pull VectorOperations into the current scope
    using Vector::abs;
    using Vector::add;
    using Vector::axpby;
    using Vector::axpy;
    using Vector::cloneVector;
    using Vector::divide;
    using Vector::dot;
    using Vector::linearSum;
    using Vector::minQuotient;
    using Vector::multiply;
    using Vector::reciprocal;
    using Vector::scale;
    using Vector::setRandomValues;
    using Vector::subtract;
    using Vector::wrmsNorm;
    using Vector::wrmsNormMask;
};


} // namespace LinearAlgebra
} // namespace AMP

#endif
