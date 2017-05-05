#ifndef included_AMP_NativeThyraVector
#define included_AMP_NativeThyraVector

#include "vectors/Vector.h"
#include "vectors/VectorEngine.h"
#include "vectors/NativeVector.h"
#include "vectors/trilinos/thyra/ThyraVector.h"


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
class NativeThyraVector : public NativeVector, public ThyraVector, public VectorEngine
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
    virtual std::string type() const  override{ return "Native Thyra Vector"; }
    virtual Vector::shared_ptr getManagedVectorCopy( AMP_MPI comm ) override;
    virtual Vector::shared_ptr getManagedVectorDuplicate( AMP_MPI comm ) override;
    virtual Vector::shared_ptr cloneVector( const Variable::shared_ptr ) const override;
    virtual void copyVector( Vector::const_shared_ptr vec ) override;
    virtual void swapVectors( Vector &other ) override;
    virtual void aliasVector( Vector & ) override;
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
    virtual void linearSum( double alpha, const VectorOperations &x, double beta, const VectorOperations &y ) override;
    virtual void axpy( double alpha, const VectorOperations &x, const VectorOperations &y ) override;
    virtual void axpby( double alpha, double beta, const VectorOperations &x ) override;
    virtual void abs( const VectorOperations &x ) override;
    virtual double min( void ) const override;
    virtual double max( void ) const override;
    virtual void setRandomValues( void ) override;
    virtual double L1Norm( void ) const override;
    virtual double L2Norm( void ) const override;
    virtual double maxNorm( void ) const override;
    virtual double dot( const VectorOperations &x ) const override;
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
    virtual BufferPtr getNewBuffer() override;
    virtual bool sameEngine( VectorEngine & ) const override;
    virtual VectorEngine::shared_ptr cloneEngine( BufferPtr p ) const override;
    virtual void swapEngines( VectorEngine::shared_ptr ) override;
    virtual void *getDataBlock( size_t i ) override;
    virtual const void *getDataBlock( size_t i ) const override;
    virtual AMP_MPI getComm() const override;
    virtual void copyOutRawData( double *out ) const override;
    virtual uint64_t getDataID() const override  { return reinterpret_cast<uint64_t>( getRawDataBlockAsVoid( 0 ) ); }


protected:
    //! Empty constructor.
    NativeThyraVector();

    virtual void *getRawDataBlockAsVoid( size_t i ) override;
    virtual const void *getRawDataBlockAsVoid( size_t i ) const override;


private:
    size_t d_local;

    static Teuchos::RCP<const Thyra::VectorBase<double>> getThyraVec( const VectorOperations &v );
    static Teuchos::RCP<const Thyra::VectorBase<double>>
    getThyraVec( const Vector::const_shared_ptr &v );


public: // Pull VectorOperations into the current scope
    using Vector::add;
    using Vector::abs;
    using Vector::axpy;
    using Vector::axpby;
    using Vector::divide;
    using Vector::dot;
    using Vector::linearSum;
    using Vector::minQuotient;
    using Vector::multiply;
    using Vector::scale;
    using Vector::setRandomValues;
    using Vector::subtract;
    using Vector::reciprocal;
    using Vector::wrmsNorm;
    using Vector::wrmsNormMask;
    using Vector::cloneVector;
};


} // LinearAlgebra namespace
} // AMP namespace


#include "NativeThyraVector.inline.h"

#endif
