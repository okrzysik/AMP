#ifndef included_AMP_EpetraVectorEngine
#define included_AMP_EpetraVectorEngine

#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include "vectors/VectorEngine.h"
#include "EpetraVector.h"


namespace AMP {
namespace LinearAlgebra {

/** \class EpetraVectorEngineParameters
  * \brief Class that details how to construct an EpetraVectorEngine
  */
class EpetraVectorEngineParameters : public VectorEngineParameters
{
public:
    /** \brief Constructor
        \param[in] local_size     The number of elements on this core
        \param[in] global_size    The number of elements in total
        \param[in] comm         Communicator to create the vector on
        \details  This assumes a contiguous allocation of data.  Core 0 has global ids \f$(0,1,\ldots,n-1)\f$, core 1 has global ids \f$(n,n+1,n+2,\ldots,m)\f$, etc.
        */
    EpetraVectorEngineParameters( size_t local_size , size_t global_size , AMP_MPI comm );

    /** \brief Constructor
      * \param[in]  local_size    The number of elements on this core
      * \param[in]  global_size   The number of elements in total
      * \param[in]  emap        An Epetra_Map for the data
      * \param[in]  ecomm       An Epetra_MpiComm for constructing the vector on
      * \details  This allows construction of an EpetraVectorEngine from handy Epetra objects
      */
    EpetraVectorEngineParameters( size_t local_size , size_t global_size , AMP::shared_ptr<Epetra_Map> emap , AMP_MPI ecomm );

    //! destructor
    virtual ~EpetraVectorEngineParameters();

    /** \brief  Return the Epetra_Map for this engine
      * \return  The Epetra_Map
      */
    Epetra_Map& getEpetraMap();

private:
    AMP::shared_ptr<Epetra_Map> d_emap;      // Epetra map
};


/** \class EpetraVectorEngine
  * \brief A linear algebra engine that uses Epetra
  * \details  Use the Epetra implementation of the L1 BLAS routines.  Unlike other
  * libraries, it is very difficult to separate the data from the engine.  For this
  * reason, the EpetraVectorEngine contains the Epetra_Vector to operate on.
  */
class EpetraVectorEngine : public VectorEngine
{
    protected:
    /** \brief  The Epetra_Vector to perform work on
      */
    Epetra_Vector  d_epetraVector;

    /** \brief  The number of local elements in the vector
      */
    int   d_iLocalSize;

    /** \brief  The number of elements in the entire vector
      */
    int   d_iGlobalSize;

    public:
    /** \brief Constructor
      * \param[in]  alias  The parameters to construct this engine
      * \param[in]  p  The buffer to use to construct the engine
      */
    EpetraVectorEngine ( VectorEngineParameters::shared_ptr  alias , BufferPtr p = BufferPtr () );
    
    /** \brief Destructor
      */
    virtual ~EpetraVectorEngine ();

    /** \brief  Get the raw Epetra_Vector
      * \return  The Epetra_Vector currently used by this engine
      */
        Epetra_Vector  &getEpetra_Vector();

    /** \brief  Get the raw Epetra_Vector
      * \return  The Epetra_Vector currently used by this engine
      */
    const Epetra_Vector  &getEpetra_Vector() const;

    AMP_MPI  getComm () const;
    virtual   void  *getDataBlock ( size_t i );
    virtual   const void  *getDataBlock ( size_t i ) const;
    virtual size_t   getLocalSize() const;
    virtual size_t   getGlobalSize() const;
    virtual BufferPtr  getNewBuffer ();
    virtual bool sameEngine ( VectorEngine &e ) const;
    virtual shared_ptr  cloneEngine ( BufferPtr p ) const;
    virtual void swapEngines ( shared_ptr );
    virtual   size_t    numberOfDataBlocks () const;
    virtual   size_t    sizeOfDataBlock ( size_t i ) const;

//    virtual void addScalar ( const VectorOperations & , double );
    virtual void setToScalar(double alpha);
    virtual void scale(double alpha, const VectorOperations &x);
    virtual void scale(double alpha);
    virtual void add(const VectorOperations &x, const VectorOperations &y);
    virtual void subtract(const VectorOperations &x, const VectorOperations &y);
    virtual void multiply(const VectorOperations &x, const VectorOperations &y);
    virtual void divide(const VectorOperations &x, const VectorOperations &y);
    virtual void reciprocal(const VectorOperations &x);
    virtual void linearSum(double alpha, const VectorOperations &x,
        double beta, const VectorOperations &y);
    virtual void axpy(double alpha, const VectorOperations &x, const VectorOperations &y);
    virtual void axpby(double alpha, double beta, const VectorOperations &x);
    virtual void abs(const VectorOperations &x);
    virtual double min(void) const;
    virtual double max(void) const;
    virtual void setRandomValues(void);
    virtual void setValuesByLocalID(int i, size_t * , const double *val);
    virtual void setLocalValuesByGlobalID(int i, size_t * , const double *val);
    virtual void addValuesByLocalID(int i, size_t * , const double *val);
    virtual void addLocalValuesByGlobalID(int i, size_t * , const double *val);
    virtual void getValuesByLocalID(int i, size_t * ,double *val) const;
    virtual void getLocalValuesByGlobalID(int i, size_t * ,double *val) const;
    double L1Norm(void) const;
    double L2Norm(void) const;
    double maxNorm(void) const;
    double dot(const VectorOperations &x) const;
    void   putRawData ( const double *in );
    void   copyOutRawData ( double *out ) const;
  };

}
}

#include "EpetraVectorEngine.inline.h"

#endif
