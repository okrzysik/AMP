#ifndef included_AMP_ManagedVector
#define included_AMP_ManagedVector


#include "AMP/vectors/DataChangeFirer.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"

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
    std::shared_ptr<VectorOperations> d_Engine;

    //! Buffer to use for the managed vector
    std::shared_ptr<VectorData> d_Buffer;
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
    virtual ~ManagedVector() {}

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
    inline std::shared_ptr<VectorOperations> getVectorEngine() { return d_Engine; }
    inline std::shared_ptr<const VectorOperations> getVectorEngine() const { return d_Engine; }

    virtual bool isAnAliasOf( Vector &rhs );
    virtual bool isAnAliasOf( Vector::shared_ptr rhs );

    virtual std::shared_ptr<ManagedVectorParameters> getManagedVectorParameters();


protected:
    //! The buffer used to store data
    std::shared_ptr<VectorData> d_vBuffer;

    //! The engine to act on the buffer
    std::shared_ptr<VectorOperations> d_Engine;

    //! The parameters used to create this vector
    std::shared_ptr<ManagedVectorParameters> d_pParameters;

    //! Function that returns a pointer to a managed vector
    virtual ManagedVector *getNewRawPtr() const = 0;


public: // Derived from VectorData
    size_t numberOfDataBlocks() const override;
    size_t sizeOfDataBlock( size_t i ) const override;
    size_t getLocalSize() const override;
    size_t getGlobalSize() const override;
    void getValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const override;
    void getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const override;
    void getGhostValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const override;
    void setValuesByGlobalID( int i, size_t *, const double *val ) override;
    void setLocalValuesByGlobalID( int i, size_t *, const double *val ) override;
    void setGhostValuesByGlobalID( int i, size_t *, const double *val ) override;
    void setValuesByLocalID( int i, size_t *, const double *val ) override;
    void addValuesByLocalID( int i, size_t *, const double *val ) override;
    void addLocalValuesByGlobalID( int i, size_t *, const double *val ) override;
    void putRawData( const double *in ) override;
    void copyOutRawData( double *in ) const override;
    UpdateState getUpdateStatus() const override;
    void setUpdateStatus( UpdateState state ) override;
    uint64_t getDataID() const override
    {
        return reinterpret_cast<uint64_t>( getRawDataBlockAsVoid( 0 ) );
    }
    bool isTypeId( size_t hash, size_t ) const override
    {
        return hash == typeid( double ).hash_code();
    }
    size_t sizeofDataBlockType( size_t ) const override { return sizeof( double ); }
    std::string VectorDataName() const override { return type(); }
    void swapData( VectorData & ) override { AMP_ERROR( "Not finished" ); }

protected: // Derived from VectorData
    void dataChanged() override;
    void *getRawDataBlockAsVoid( size_t i ) override;
    const void *getRawDataBlockAsVoid( size_t i ) const override;


public: // Derived from VectorOperations
    double L1Norm( void ) const override;
    double L2Norm( void ) const override;
    double maxNorm( void ) const override;
    double dot( const VectorOperations &x ) const override;
    void axpy( double alpha, const VectorOperations &x, const VectorOperations &y ) override;
    void axpby( double alpha, double beta, const VectorOperations &x ) override;
    void abs( const VectorOperations &x ) override;
    double min( void ) const override;
    double max( void ) const override;
    void setRandomValues( void ) override;
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


public: // Derived from Vector
    using Vector::cloneVector;
    std::string type() const override;
    std::shared_ptr<Vector> cloneVector( const Variable::shared_ptr name ) const override;
    std::shared_ptr<ParameterBase> getParameters() override;
    Vector::shared_ptr subsetVectorForVariable( Variable::const_shared_ptr name ) override;
    Vector::const_shared_ptr
    constSubsetVectorForVariable( Variable::const_shared_ptr name ) const override;
    void copy( const VectorOperations &src ) override;
    void swapVectors( Vector &other ) override;
    void aliasVector( Vector &other ) override;


protected: // Derived from Vector
    Vector::shared_ptr selectInto( const VectorSelector & ) override;
    Vector::const_shared_ptr selectInto( const VectorSelector & ) const override;
    void addCommunicationListToParameters( CommunicationList::shared_ptr comm ) override;


public: // Pull VectorOperations into the current scope
    using VectorOperations::abs;
    using VectorOperations::add;
    using VectorOperations::addScalar;
    using VectorOperations::axpby;
    using VectorOperations::axpy;
    using VectorOperations::divide;
    using VectorOperations::dot;
    using VectorOperations::equals;
    using VectorOperations::linearSum;
    using VectorOperations::minQuotient;
    using VectorOperations::multiply;
    using VectorOperations::reciprocal;
    using VectorOperations::scale;
    using VectorOperations::setRandomValues;
    using VectorOperations::subtract;
    using VectorOperations::wrmsNorm;
    using VectorOperations::wrmsNormMask;
    using VectorOperations::zero;

private:
    ManagedVector();
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
