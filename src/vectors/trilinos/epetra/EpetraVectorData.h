#ifndef included_AMP_EpetraVectorData
#define included_AMP_EpetraVectorData

#include "vectors/data/VectorData.h"

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

    virtual std::string VectorDataName() const override { return "EpetraVectorData"; }
    virtual size_t getLocalSize() const override { return d_iLocalSize; }
    virtual size_t getGlobalSize() const override { return d_iGlobalSize; }
    virtual size_t getLocalStartID() const override { return d_iLocalStart; }
    virtual size_t numberOfDataBlocks() const override { return 1; }
    virtual size_t sizeOfDataBlock( size_t i ) const override { return i == 0 ? d_iLocalSize : 0; }
    virtual void setValuesByLocalID( int i, size_t *, const double *val ) override;
    virtual void setLocalValuesByGlobalID( int i, size_t *, const double *val ) override;
    virtual void addValuesByLocalID( int i, size_t *, const double *val ) override;
    virtual void addLocalValuesByGlobalID( int i, size_t *, const double *val ) override;
    virtual void getValuesByLocalID( int i, size_t *, double *val ) const override;
    virtual void getLocalValuesByGlobalID( int i, size_t *, double *val ) const override;
    virtual void putRawData( const double *in ) override;
    virtual void copyOutRawData( double *out ) const override;
    virtual uint64_t getDataID() const override
    {
        return reinterpret_cast<uint64_t>( getRawDataBlockAsVoid( 0 ) );
    }
    virtual void *getRawDataBlockAsVoid( size_t i ) override;
    virtual const void *getRawDataBlockAsVoid( size_t i ) const override;
    virtual size_t sizeofDataBlockType( size_t ) const override { return sizeof( double ); }
    virtual bool isTypeId( size_t hash, size_t ) const override
    {
        return hash == typeid( double ).hash_code();
    }

protected:
    EpetraVectorData( Epetra_DataAccess, const Epetra_BlockMap &, double *, int, int, int );

    //! The Epetra_Vector to perform work on
    Epetra_Vector d_epetraVector;

    //! A shared ptr to the buffer (to prevent it from leaving scope)
    AMP::shared_ptr<VectorData> d_buf_scope;

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
