#ifndef included_AMP_ManagedVectorData
#define included_AMP_ManagedVectorData

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/DataChangeListener.h"
#include "AMP/vectors/data/VectorData.h"

#include <stdexcept>
#include <vector>


namespace AMP::LinearAlgebra {


/**
   \brief Class used to control data and kernels of various vector libraries
   \details  A ManagedVector will take an engine and create a buffer, if
   necessary.

   A ManagedVector has two pointers: data and engine.  If the data pointer
   is null, then the engine is assumed to have the data.
*/
class ManagedVectorData : public VectorData, public DataChangeListener
{

public:
    /** \brief Construct a ManagedVector from a set of parameters
     * \param[in] params  The description of the ManagedVector
     */
    explicit ManagedVectorData( std::shared_ptr<Vector> vec );

    /** \brief Construct a view of an AMP vector
     * \param[in] alias  Vector to view
     */
    explicit ManagedVectorData( std::shared_ptr<VectorData> alias );

    //! Destructor
    virtual ~ManagedVectorData();

    /** \brief  Return the engine associated with this ManagedVector
     * \return The engine
     */
    std::shared_ptr<Vector> getVectorEngine();
    std::shared_ptr<const Vector> getVectorEngine() const;

    bool isAnAliasOf( const VectorData &rhs ) const override;

    void receiveDataChanged() override { fireDataChange(); }

protected:
    //! The unterlying vector
    std::shared_ptr<Vector> d_Engine;


public: // Derived from VectorData
    size_t numberOfDataBlocks() const override;
    size_t sizeOfDataBlock( size_t i ) const override;
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
    std::string VectorDataName() const override;
    void swapData( VectorData & ) override;

    void dataChanged() override;

    std::shared_ptr<VectorData> cloneData() const override;

protected: // Derived from VectorData
    void *getRawDataBlockAsVoid( size_t i ) override;
    const void *getRawDataBlockAsVoid( size_t i ) const override;

private:
    ManagedVectorData();
};


} // namespace AMP::LinearAlgebra


#endif
