#ifndef  included_AMP_SimpleVector
#define  included_AMP_SimpleVector

#include "Vector.h"

namespace AMP {
namespace LinearAlgebra {


/** \brief A core-local vector
* \details This is a Vector that implements the Vector interface for a std::vector<double>.
*/
template <typename T>
class SimpleVector : public Vector
{
protected:
    std::vector<T>  d_Data;
    size_t d_startIndex;
    size_t d_globalSize;
    AMP_MPI d_comm;

    SimpleVector();
    SimpleVector( const SimpleVector & );

public:
    /** \brief    Create a SimpleVector
      * \details  This is the factory method for the SimpleVector.  It returns the shared pointer
      * to be used in the code
      * \param    localSize  The number of elements in the vector on this processor
      * \param    var The variable associated with the new vector
      */
    static Vector::shared_ptr  create( size_t localSize, Variable::shared_ptr var );

    /** \brief    Create a SimpleVector
      * \details  This is the factory method for the SimpleVector.  It returns the shared pointer
      * to be used in the code
      * \param    localSize  The number of elements in the vector on this processor
      * \param    var The variable associated with the new vector
      * \param    comm The variable associated with the new vector
      */
    static Vector::shared_ptr  create( size_t localSize, Variable::shared_ptr var, AMP_MPI comm );

    /** \brief    Create a SimpleVector
      * \details  This is the factory method for the SimpleVector.  It returns the shared pointer
      * to be used in the code that spans a comm and contains ghost values.
      * \param    var The variable associated with the new vector
      * \param    DOFs The DOFManager
      * \param    commlist The communication list
      */
    static Vector::shared_ptr  create(  Variable::shared_ptr var,
      AMP::Discretization::DOFManager::shared_ptr DOFs, 
      AMP::LinearAlgebra::CommunicationList::shared_ptr commlist );

    /** \brief  Destructor
      */
    virtual ~SimpleVector() override {}

    virtual std::string type() const override { return "Simple Vector"; }
    virtual uint64_t getDataID() const override { return reinterpret_cast<uint64_t>(d_Data.data()); }
    virtual void setToScalar(double alpha) override;
    virtual void scale(double alpha, const VectorOperations &x) override;
    virtual void scale(double alpha) override;
    virtual void add(const VectorOperations &x, const VectorOperations &y) override;
    virtual void subtract(const VectorOperations &x, const VectorOperations &y) override;
    virtual void multiply( const VectorOperations &x, const VectorOperations &y) override;
    virtual void divide( const VectorOperations &x, const VectorOperations &y) override;
    virtual void reciprocal(const VectorOperations &x) override;
    virtual void linearSum(double alpha, const VectorOperations &x,
          double beta, const VectorOperations &y) override;
    virtual void axpy(double alpha, const VectorOperations &x, const VectorOperations &y) override;
    virtual void axpby(double alpha, double beta, const VectorOperations &x) override;
    virtual void abs(const VectorOperations &x) override;
    virtual double min(void) const override;
    virtual double max(void) const override;
    virtual double L1Norm(void) const override;
    virtual double L2Norm(void) const override;
    virtual double maxNorm(void) const override;

    using Vector::dot;
    virtual double dot(const VectorOperations &x) const override;

    virtual AMP::shared_ptr<ParameterBase> getParameters ();
    using Vector::cloneVector;
    virtual Vector::shared_ptr cloneVector(const Variable::shared_ptr name) const override;
    virtual size_t  numberOfDataBlocks() const override;
    virtual size_t  sizeOfDataBlock( size_t i = 0 ) const override;
    using Vector::copyVector;
    virtual void copyVector( Vector::const_shared_ptr src_vec ) override;
    virtual void swapVectors(Vector &other) override;
    virtual void aliasVector(Vector &other) override;
    virtual void setValuesByLocalID( int num , size_t *indices , const double *vals ) override;

    /** \brief Not implemented
      */
    virtual void setLocalValuesByGlobalID( int num , size_t *indices , const double *vals ) override;
    virtual void addValuesByLocalID( int num , size_t *indices , const double *vals ) override;

    /** \brief Not implemented
      */
    virtual void addLocalValuesByGlobalID( int num , size_t *indices , const double *vals ) override;

    /** \brief Not implemented
      */
    virtual void getLocalValuesByGlobalID( int num , size_t *indices , double *vals ) const override;
    virtual void assemble() override;
    virtual void putRawData( const double *in ) override;
    virtual void copyOutRawData( double *out ) const override;
    virtual size_t getLocalSize() const override;
    virtual size_t getGlobalSize() const override;
    virtual void *getRawDataBlockAsVoid( size_t i ) override;
    virtual const void *getRawDataBlockAsVoid( size_t i ) const override;

    T &operator[] ( size_t i );
    T operator[] ( size_t i ) const ;

    /** \brief Resize this vector
      * \param[in] i The new size
      */
    void  resize( size_t i );

    //! return a const reference to the internal data container
    const std::vector<T> &getData( void ) const { return d_Data; }
};


}
}

#include "SimpleVector.inline.h"
#include "SimpleVector.hpp"

#endif
