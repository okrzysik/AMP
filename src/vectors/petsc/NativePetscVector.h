#ifndef included_AMP_NativePetscVector
#define included_AMP_NativePetscVector

#include "vectors/NativeVector.h"
#include "vectors/petsc/ManagedPetscVector.h"
#include "utils/AMP_MPI.h"

extern "C"{
#include "petsc.h"
#include "petscvec.h"
}


namespace AMP {
namespace LinearAlgebra {



/** \class NativePetscVectorParameters
  * \brief Parameters to set when creating a NativePetscVector
  */
class NativePetscVectorParameters : public NativeVectorParameters
{
public:
    //!  The vector to wrap
    Vec         d_InVec;

    //!  The communicator associated with the Vec
    AMP_MPI     d_Comm;

    //!  If true, ~NativePetscVector() will call VecDestroy()
    bool        d_Deleteable;

    //! The number of local entities in the vector
    size_t      d_localsize;

    /** \brief Constructor
      * \param[in] v The vector to wrap
      */
    NativePetscVectorParameters ( Vec v );
};



/** \class NativePetscVector
  * \brief An AMP Vector that uses PETSc for parallel data management, linear algebra,
  * etc.
  * \details  This is an AMP wrapper to PETSc.  This is different from ManagedPetscVector
  * in that this class does not replace calls to Vec*.  Rather, it wraps these calls.
  * This class is used when PETSc is chosen as the default linear algebra engine.
  *
  * This class is not to be used directly, just through base class interfaces.
  * \see PetscVector
  * \see ManagedPetscVector
  */
class NativePetscVector : public NativeVector , 
                          public PetscVector , 
                          public VectorEngine
{
public:
      /** \brief Conveninece typedef
        */
      typedef  NativeVector::parameters_ptr           parameters_ptr;

      /** \brief Conveninece typedef
        */
      typedef  NativeVectorParameters                 parameters;


      /** \brief Construct a wrapper for a PETSc Vec from a set of parameters
        * \param[in] params The parameters describing the Vec
        */
      NativePetscVector ( VectorParameters::shared_ptr params );
      /** \brief Destructor
        */
      virtual ~NativePetscVector ();

      virtual std::string type() const { return "Native PETSc Vector"; }

      virtual Vector::shared_ptr  getManagedVectorCopy ( AMP_MPI comm );
 
      virtual Vector::shared_ptr  getManagedVectorDuplicate ( AMP_MPI  comm );

      virtual Vector::shared_ptr cloneVector(const Variable::shared_ptr ) const;
      virtual void copyVector(const Vector &src_vec);

      virtual void swapVectors(Vector &other);
      virtual void aliasVector(Vector & );

      virtual   size_t      numberOfDataBlocks () const;
      virtual   size_t      sizeOfDataBlock ( size_t i ) const;


      virtual void setToScalar(double alpha);
      virtual void scale(double alpha, const VectorOperations &x);
      virtual void scale(double alpha);
      virtual void add(const VectorOperations &x, const VectorOperations &y);
      virtual void subtract(const VectorOperations &x, const VectorOperations &y);
      virtual void multiply( const VectorOperations &x, const VectorOperations &y);
      virtual void divide( const VectorOperations &x, const VectorOperations &y);
      virtual void reciprocal(const VectorOperations &x);
      virtual void linearSum(double alpha, const VectorOperations &x,
              double beta, const VectorOperations &y);
      virtual void axpy(double alpha, const VectorOperations &x, const VectorOperations &y);
      virtual void axpby(double alpha, double beta, const VectorOperations &x);
      virtual void abs(const VectorOperations &x);
      virtual double min(void) const;
      virtual double max(void) const;
      virtual void setRandomValues(void);

      virtual double L1Norm(void) const;
      virtual double L2Norm(void) const;
      virtual double maxNorm(void) const;
      virtual double dot(const VectorOperations &x) const;
    
      virtual double localL1Norm(void) const;
      virtual double localL2Norm(void) const;
      virtual double localMaxNorm(void) const;

      virtual void setValuesByLocalID(int , size_t * , const double *);
      virtual void setLocalValuesByGlobalID(int , size_t * , const double *);
      virtual void addValuesByLocalID(int , size_t * , const double *);
      virtual void addLocalValuesByGlobalID(int , size_t * , const double *);

      virtual void getLocalValuesByGlobalID ( int numVals , size_t *ndx , double *vals ) const;
      virtual void getValuesByLocalID ( int numVals , size_t *ndx , double *vals ) const;

      virtual void assemble();

      virtual size_t getLocalSize() const;
      virtual size_t getGlobalSize() const;

      virtual void putRawData ( double * );

      virtual   BufferPtr   getNewBuffer ();
      virtual   bool        sameEngine ( VectorEngine & ) const;
      virtual   VectorEngine::shared_ptr  cloneEngine ( BufferPtr p ) const;

      virtual   void        swapEngines ( VectorEngine::shared_ptr );

      virtual   void       *getDataBlock ( size_t i );
      virtual   const void *getDataBlock ( size_t i ) const;

      virtual   AMP_MPI   getComm() const;

      virtual void copyOutRawData ( double **out );

      virtual boost::shared_ptr<ParameterBase> getParameters ();

protected:

      void *getRawDataBlockAsVoid ( size_t i );
      const void *getRawDataBlockAsVoid ( size_t i ) const;

      void  resetArray ();
      void  resetArray () const;

private:
      parameters_ptr   d_pParameters;
      bool              d_bDeleteMe;
      mutable double   *d_pArray;  //mutable so that we can cache the value
};


}
}

#include "NativePetscVector.inline.h"

#endif
