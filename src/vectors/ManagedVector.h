#ifndef included_AMP_ManagedVector
#define included_AMP_ManagedVector


#include "vectors/DataChangeFirer.h"
#include "vectors/Vector.h"
#include "vectors/VectorEngine.h"
#include "vectors/operations/VectorOperationsDefault.h"

#include <stdexcept>
#include <vector>


namespace AMP {
namespace LinearAlgebra {


/**
  \brief Data necessary to create a managed vector
*/
class ManagedVectorParameters : public VectorParameters
{
protected:
    //!  Copy constructor is protected to prevent unintended copies
    ManagedVectorParameters( const ManagedVectorParameters & );

public:
    //! Constructor
    ManagedVectorParameters();

    //! The VectorEngine to use with the managed vector
    VectorEngine::shared_ptr d_Engine;

    //! Indicates whether the engine should be used as is or cloned
    bool d_CloneEngine;

    //! Buffer to use for the managed vector
    VectorEngine::BufferPtr d_Buffer;
};


/**
   \brief Class used to control data and kernels of various vector libraries
   \details  A ManagedVector will take an engine and create a buffer, if
   necessary.

   A ManagedVector has two pointers: data and engine.  If the data pointer
   is null, then the engine is assumed to have the data.
*/
class ManagedVector : public Vector, public VectorOperationsDefault<double>, public DataChangeFirer
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
    VectorEngine::shared_ptr getVectorEngine();
    VectorEngine::const_shared_ptr getVectorEngine() const;
    std::string type() const override;

    virtual Vector::shared_ptr subsetVectorForVariable( Variable::const_shared_ptr name ) override;
    virtual Vector::const_shared_ptr
    constSubsetVectorForVariable( Variable::const_shared_ptr name ) const override;
    virtual size_t numberOfDataBlocks() const override;
    virtual size_t sizeOfDataBlock( size_t i ) const override;
    virtual void copy( const VectorOperations &src ) override;
    virtual void swapVectors( Vector &other ) override;
    virtual void aliasVector( Vector &other ) override;

    virtual bool isAnAliasOf( Vector &rhs );
    virtual bool isAnAliasOf( Vector::shared_ptr rhs );
    virtual AMP::shared_ptr<Vector> cloneVector( const Variable::shared_ptr name ) const override;
    virtual AMP::shared_ptr<ParameterBase> getParameters() override;

    virtual AMP::shared_ptr<ManagedVectorParameters> getManagedVectorParameters();

    virtual size_t getLocalSize() const override;
    virtual size_t getGlobalSize() const override;

    virtual void getValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const override;
    virtual void getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const override;
    virtual void getGhostValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const override;
    virtual void setValuesByGlobalID( int i, size_t *, const double *val ) override;
    virtual void setLocalValuesByGlobalID( int i, size_t *, const double *val ) override;
    virtual void setGhostValuesByGlobalID( int i, size_t *, const double *val ) override;

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
    virtual void setValuesByLocalID( int i, size_t *, const double *val ) override;
    virtual void addValuesByLocalID( int i, size_t *, const double *val ) override;
    virtual void addLocalValuesByGlobalID( int i, size_t *, const double *val ) override;
    virtual void putRawData( const double *in ) override;
    virtual void copyOutRawData( double *in ) const override;
    double L1Norm( void ) const override;
    double L2Norm( void ) const override;
    double maxNorm( void ) const override;
    double dot( const VectorOperations &x ) const override;
    virtual UpdateState getUpdateStatus() const override;
    virtual void setUpdateStatus( UpdateState state ) override;
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
    virtual Vector::shared_ptr selectInto( const VectorSelector & ) override;
    virtual Vector::const_shared_ptr selectInto( const VectorSelector & ) const override;


    //! The buffer used to store data
    VectorEngine::BufferPtr d_vBuffer;

    //! The engine to act on the buffer
    VectorEngine::shared_ptr d_Engine;

    //! The parameters used to create this vector
    AMP::shared_ptr<ManagedVectorParameters> d_pParameters;

    //! Function that returns a pointer to a managed vector
    virtual ManagedVector *getNewRawPtr() const = 0;

    virtual void dataChanged() override;
    virtual void *getRawDataBlockAsVoid( size_t i ) override;
    virtual const void *getRawDataBlockAsVoid( size_t i ) const override;

    virtual void addCommunicationListToParameters( CommunicationList::shared_ptr comm ) override;


public: // Pull Vector into the current scope
    using Vector::cloneVector;

public: // Pull VectorOperations into the current scope
    using VectorOperations::add;
    using VectorOperations::addScalar;
    using VectorOperations::abs;
    using VectorOperations::axpy;
    using VectorOperations::axpby;
    using VectorOperations::divide;
    using VectorOperations::dot;
    using VectorOperations::equals;
    using VectorOperations::linearSum;
    using VectorOperations::minQuotient;
    using VectorOperations::multiply;
    using VectorOperations::setRandomValues;
    using VectorOperations::scale;
    using VectorOperations::subtract;
    using VectorOperations::reciprocal;
    using VectorOperations::wrmsNorm;
    using VectorOperations::wrmsNormMask;
    using VectorOperations::zero;

private:
    ManagedVector();
};
}
}

#include "ManagedVector.inline.h"

#endif
