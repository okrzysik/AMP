#ifndef included_AMP_NativeThyraVector
#define included_AMP_NativeThyraVector

#include "vectors/NativeVector.h"
#include "vectors/trilinos/ThyraVector.h"


namespace AMP {
namespace LinearAlgebra {



/** \class NativeThyraVectorParameters
  * \brief Parameters to set when creating a NativeThyraVector
  */
class NativeThyraVectorParameters : public NativeVectorParameters
{
public:

    //! The vector to wrap
    Teuchos::RCP<Thyra::VectorSpaceBase<double> >  d_InVec;

    /** \brief Constructor
      * \details  This will create NativeThyraVectorParameters that can be used
      *   to construct a NativeThyraVector around the given vector.
      * \param[in] v            The vector to wrap
      */
    NativeThyraVectorParameters( Teuchos::RCP<Thyra::VectorSpaceBase<double> > v );
};



/** \class NativeThyraVector
  * \brief An AMP Vector that uses Thyra for parallel data management, linear algebra,
  * etc.
  * \details  This is an AMP wrapper to Thyra.  This is different from ManagedThyraVector
  * in that this class does not replace calls to Vec*.  Rather, it wraps these calls.
  * This class is used when Thyra is chosen as the default linear algebra engine.
  *
  * This class is not to be used directly, just through base class interfaces.
  * \see ThyraVector
  * \see ManagedThyraVector
  */
class NativeThyraVector : public NativeVector , 
                          public ThyraVector , 
                          public VectorEngine
{
public:
      /** \brief Conveninece typedef
        */
      typedef  NativeVector::parameters_ptr           parameters_ptr;

      /** \brief Conveninece typedef
        */
      typedef  NativeVectorParameters                 parameters;


      /** \brief Construct a wrapper for a Thyra Vec from a set of parameters
        * \param[in] params The parameters describing the Vec
        */
      NativeThyraVector ( VectorParameters::shared_ptr params );
      /** \brief Destructor
        */
      virtual ~NativeThyraVector ();

      virtual std::string type() const { return "Native Thyra Vector"; }

      virtual Vector::shared_ptr  getManagedVectorCopy ( AMP_MPI comm );
 
      virtual Vector::shared_ptr  getManagedVectorDuplicate ( AMP_MPI  comm );

      using Vector::cloneVector;
      virtual Vector::shared_ptr cloneVector(const Variable::shared_ptr ) const;
      virtual void copyVector(const Vector::const_shared_ptr &vec);

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
      using Vector::dot;
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

protected:

private:

};


}
}


#include "NativeThyraVector.inline.h"

#endif
