#ifndef included_AMP_SimpleVector
#define included_AMP_SimpleVector

#include "Vector.h"
#include "vectors/operations/VectorOperationsDefault.h"
#include "vectors/operations/VectorOperationsDefault.hpp"


namespace AMP {
namespace LinearAlgebra {


/** \brief A core-local vector
* \details This is a Vector that implements the Vector interface for a std::vector<double>.
*/
template <typename T>
class SimpleVector :
    public Vector,
    public VectorOperationsDefault<T>
{
protected:
    std::vector<T> d_Data;
    size_t d_startIndex;
    size_t d_globalSize;
    AMP_MPI d_comm;

    SimpleVector();
    explicit SimpleVector( const SimpleVector & );

public:
    /** \brief    Create a SimpleVector
      * \details  This is the factory method for the SimpleVector.  It returns the shared pointer
      * to be used in the code
      * \param    localSize  The number of elements in the vector on this processor
      * \param    var The variable associated with the new vector
      */
    static Vector::shared_ptr create( size_t localSize, const std::string& var );

    /** \brief    Create a SimpleVector
      * \details  This is the factory method for the SimpleVector.  It returns the shared pointer
      * to be used in the code
      * \param    localSize  The number of elements in the vector on this processor
      * \param    var The variable associated with the new vector
      */
    static Vector::shared_ptr create( size_t localSize, Variable::shared_ptr var );

    /** \brief    Create a SimpleVector
      * \details  This is the factory method for the SimpleVector.  It returns the shared pointer
      * to be used in the code
      * \param    localSize  The number of elements in the vector on this processor
      * \param    var The variable associated with the new vector
      * \param    comm The variable associated with the new vector
      */
    static Vector::shared_ptr create( size_t localSize, Variable::shared_ptr var, AMP_MPI comm );

    /** \brief    Create a SimpleVector
      * \details  This is the factory method for the SimpleVector.  It returns the shared pointer
      * to be used in the code that spans a comm and contains ghost values.
      * \param    var The variable associated with the new vector
      * \param    DOFs The DOFManager
      * \param    commlist The communication list
      */
    static Vector::shared_ptr create( Variable::shared_ptr var,
                                      AMP::Discretization::DOFManager::shared_ptr DOFs,
                                      AMP::LinearAlgebra::CommunicationList::shared_ptr commlist );

    /** \brief  Destructor
      */
    virtual ~SimpleVector() override {}

    virtual std::string type() const override { return "Simple Vector"; }
    virtual uint64_t getDataID() const override
    {
        return reinterpret_cast<uint64_t>( d_Data.data() );
    }
    virtual Vector::shared_ptr cloneVector( const Variable::shared_ptr name ) const override;
    virtual size_t numberOfDataBlocks() const override;
    virtual size_t sizeOfDataBlock( size_t i = 0 ) const override;
    virtual void copyVector( Vector::const_shared_ptr src_vec ) override;
    virtual void swapVectors( Vector &other ) override;
    virtual void aliasVector( Vector &other ) override;
    virtual void setValuesByLocalID( int num, size_t *indices, const double *vals ) override;

    /** \brief Not implemented
      */
    virtual void setLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) override;
    virtual void addValuesByLocalID( int num, size_t *indices, const double *vals ) override;

    /** \brief Not implemented
      */
    virtual void addLocalValuesByGlobalID( int num, size_t *indices, const double *vals ) override;

    /** \brief Not implemented
      */
    virtual void getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const override;
    virtual void assemble() override;
    virtual void putRawData( const double *in ) override;
    virtual void copyOutRawData( double *out ) const override;
    virtual size_t getLocalSize() const override;
    virtual size_t getGlobalSize() const override;
    virtual void *getRawDataBlockAsVoid( size_t i ) override;
    virtual const void *getRawDataBlockAsVoid( size_t i ) const override;

    T &operator[]( size_t i );
    T operator[]( size_t i ) const;

    /** \brief Resize this vector
      * \param[in] i The new size
      */
    virtual void resize( size_t i );

    //! return a const reference to the internal data container
    const std::vector<T> &getData( void ) const { return d_Data; }

protected:
    virtual bool isTypeId( size_t hash, size_t ) const override { return hash == typeid(T).hash_code(); }

public:
    using Vector::cloneVector;
    using Vector::copyVector;

};


}
}

#include "SimpleVector.hpp"
#include "SimpleVector.inline.h"

#endif
