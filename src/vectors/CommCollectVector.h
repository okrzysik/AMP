#ifndef included_AMP_CommCollectVector
#define included_AMP_CommCollectVector

#include "Vector.h"
#include "VectorEngine.h"
#include "DataChangePassThrough.h"
#include "utils/AMP_MPI.h"

namespace AMP {
namespace LinearAlgebra {

  typedef VectorParameters  CommCollectVectorParameters;

  class CommCollectVector : public Vector , public VectorEngine , public DataChangePassThrough
  {
    protected:
      Vector::shared_ptr    d_SmallCommVector;
      AMP_MPI               d_LargeComm;
      AMP_MPI               d_RankComm;

      size_t                d_LocalSize;
      size_t                d_GlobalSize;

      void  makeRankComm ();
      VectorEngine::shared_ptr  getVectorEngine();
      const VectorEngine::shared_ptr  getVectorEngine() const;

      virtual void *getRawDataBlockAsVoid ( size_t i );

      virtual const void *getRawDataBlockAsVoid ( size_t i ) const;
      const VectorOperations &getOps ( const VectorOperations &ops ) const;
      VectorOperations &getOps ( VectorOperations &ops ) const;
      static Vector::shared_ptr  internalView ( Vector::shared_ptr p , AMP_MPI lc , AMP_MPI rc );


    public:

      Vector::shared_ptr   getSmallCommVector ();
      const Vector::shared_ptr   getSmallCommVector () const;

      static Vector::shared_ptr  view ( Vector::shared_ptr p , AMP_MPI c=AMP_MPI(AMP_COMM_WORLD) );

      static const Vector::shared_ptr  constView ( Vector::shared_ptr p , AMP_MPI c=AMP_MPI(AMP_COMM_WORLD) );


      virtual Vector::iterator    begin();
      virtual Vector::iterator    end();
      virtual Vector::const_iterator    begin() const;
      virtual Vector::const_iterator    end() const;

      virtual void        selectInto ( const VectorSelector &a , Vector::shared_ptr b );
      virtual void        dumpOwnedData ( std::ostream &out , size_t GIDoffset=0 , size_t LIDoffset = 0 ) const;
      virtual void        dumpGhostedData ( std::ostream &out , size_t offset=0 ) const;
      virtual std::string type() const;
      virtual AMP_MPI     getComm () const;
      virtual size_t numberOfDataBlocks () const;
      virtual size_t sizeOfDataBlock ( size_t i ) const;

      virtual Vector::shared_ptr  subsetVectorForVariable ( const Variable::shared_ptr &name );
      virtual Vector::shared_ptr cloneVector(const Variable::shared_ptr name) const;
      virtual void swapVectors(Vector &other);
      virtual void aliasVector(Vector &other);
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
      virtual void setRandomValues(void);
      virtual void setValuesByLocalID ( int num , int *indices , const double *vals );
      virtual void setLocalValuesByGlobalID ( int num , int *indices , const double *vals );
      virtual void setValuesByGlobalID ( int num , int *indices , const double *vals );
      virtual void addValuesByLocalID ( int num , int *indices , const double *vals );
      virtual void addLocalValuesByGlobalID ( int num , int *indices , const double *vals );
      virtual void addValuesByGlobalID ( int num , int *indices , const double *vals );
      virtual void getValuesByGlobalID ( int numVals , int *ndx , double *vals ) const;
      virtual void getLocalValuesByGlobalID ( int numVals , int *ndx , double *vals ) const;
      virtual void getValuesByLocalID ( int numVals , int *ndx , double *vals ) const;
      virtual void makeConsistent ( ScatterType  t );
      virtual void assemble();
      virtual size_t getLocalSize() const;
      virtual size_t getGlobalSize() const;
      virtual size_t getGhostSize() const;
      virtual void   putRawData ( double *d );
      virtual double min(void) const;
      virtual double max(void) const;
      virtual double L1Norm(void) const;
      virtual double L2Norm(void) const;
      virtual double maxNorm(void) const;
      virtual double dot(const VectorOperations &x) const;

      ///Vector engine functions
      virtual VectorEngine::BufferPtr  getNewBuffer();
      virtual bool                     sameEngine ( VectorEngine &rhs ) const;
      virtual VectorEngine::shared_ptr cloneEngine ( VectorEngine::BufferPtr p ) const;
      virtual void                     swapEngines ( VectorEngine::shared_ptr p );
      virtual const void              *getDataBlock ( size_t i ) const;
      virtual void                    *getDataBlock ( size_t i );
      virtual void   copyOutRawData ( double ** out );
  };

}
}

#include "CommCollectVector.inline.h"

#endif
