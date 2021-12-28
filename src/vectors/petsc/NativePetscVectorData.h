#ifndef included_AMP_NativePetscVectorData
#define included_AMP_NativePetscVectorData

#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/data/VectorData.h"

#include "petscvec.h"


namespace AMP::LinearAlgebra {


/*! \struct Vec
    \brief PETSc vector
*/


/** \class NativePetscVector
 * \brief An AMP Vector that uses PETSc for parallel data management, linear algebra,
 * etc.
 * \details  This is an AMP wrapper to PETSc.   Rather, it wraps these calls.
 * This class is used when PETSc is chosen as the default linear algebra engine.
 */
class NativePetscVectorData : public VectorData
{
public:
    /** \brief Construct a wrapper for a PETSc Vec from a set of parameters
     * \param[in] v             The Vec to wrap
     * \param[in] deleteable    If true, ~NativePetscVector() will call VecDestroy()
     * \param[in] comm          The communicator associated with the Vec
     */
    explicit NativePetscVectorData( Vec v, bool deleteable, AMP_MPI comm = AMP_MPI() );

    //! Destructor
    virtual ~NativePetscVectorData();

    std::string VectorDataName() const override { return "NativePetscVector"; }
    size_t numberOfDataBlocks() const override;
    size_t sizeOfDataBlock( size_t i ) const override;
    void putRawData( const double * ) override;
    void copyOutRawData( double *out ) const override;
    size_t getLocalSize() const override;
    size_t getGlobalSize() const override;

    void setValuesByLocalID( int, size_t *, const double * ) override;
    void setLocalValuesByGlobalID( int, size_t *, const double * ) override;
    void addValuesByLocalID( int, size_t *, const double * ) override;
    void addLocalValuesByGlobalID( int, size_t *, const double * ) override;

    void getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const override;

    // Return the id of the data
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
    void swapData( VectorData &rhs ) override;

    /** \brief Clone the data
     */
    std::shared_ptr<VectorData> cloneData() const override;

    void assemble() override;

    //! Get the PETSc vector
    Vec &getVec() { return d_petscVec; }

    //! Get the PETSc vector
    const Vec &getVec() const { return d_petscVec; }

protected:
    void resetArray();
    void resetArray() const;

private:
    friend class NativePetscVectorOperations;
    bool d_bDeleteMe;
    Vec d_petscVec;
    mutable double *d_pArray; // mutable so that we can cache the value
};


} // namespace AMP::LinearAlgebra

#endif
