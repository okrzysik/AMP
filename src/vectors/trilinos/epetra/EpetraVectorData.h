#ifndef included_AMP_EpetraVectorData
#define included_AMP_EpetraVectorData

#include "AMP/vectors/data/VectorData.h"

#include <Epetra_Vector.h>


namespace AMP {
namespace LinearAlgebra {


/** \class EpetraVectorData
 * \brief A linear algebra engine that uses Epetra
 * \details  Use the Epetra implementation of the L1 BLAS routines.  Unlike other
 * libraries, it is very difficult to separate the data from the engine.  For this
 * reason, the EpetraVectorEngine contains the Epetra_Vector to operate on.
 */
class EpetraVectorData : virtual public VectorData
{

public: // Virtual functions
    //! Virtual destructor
    virtual ~EpetraVectorData() {}

    std::string VectorDataName() const override { return "EpetraVectorData"; }
    size_t getLocalSize() const override { return d_iLocalSize; }
    size_t getGlobalSize() const override { return d_iGlobalSize; }
    size_t getLocalStartID() const override { return d_iLocalStart; }
    size_t numberOfDataBlocks() const override { return 1; }
    size_t sizeOfDataBlock( size_t i ) const override { return i == 0 ? d_iLocalSize : 0; }
    void setValuesByLocalID( int i, size_t *, const double *val ) override;
    void setLocalValuesByGlobalID( int i, size_t *, const double *val ) override;
    void addValuesByLocalID( int i, size_t *, const double *val ) override;
    void addLocalValuesByGlobalID( int i, size_t *, const double *val ) override;
    void getValuesByLocalID( int i, size_t *, double *val ) const override;
    void getLocalValuesByGlobalID( int i, size_t *, double *val ) const override;
    void putRawData( const double *in ) override;
    void copyOutRawData( double *out ) const override;
    uint64_t getDataID() const override
    {
        return reinterpret_cast<uint64_t>( getRawDataBlockAsVoid( 0 ) );
    }
    void *getRawDataBlockAsVoid( size_t i ) override;
    const void *getRawDataBlockAsVoid( size_t i ) const override;
    size_t sizeofDataBlockType( size_t ) const override { return sizeof( double ); }
    bool isTypeId( size_t hash, size_t ) const override
    {
        return hash == typeid( double ).hash_code();
    }
    void swapData( VectorData & ) override { AMP_ERROR( "Not finished" ); }
    std::shared_ptr<VectorData> cloneData() const override;


protected:
    EpetraVectorData( Epetra_DataAccess, const Epetra_BlockMap &, double *, int, int, int );

    //! The Epetra_Vector to perform work on
    Epetra_Vector d_epetraVector;

    //! A shared ptr to the buffer (to prevent it from leaving scope)
    std::shared_ptr<VectorData> d_buf_scope;

    //! The local start index
    int d_iLocalStart;

    //! The number of local elements in the vector
    int d_iLocalSize;

    //! The number of elements in the entire vector
    int d_iGlobalSize;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
