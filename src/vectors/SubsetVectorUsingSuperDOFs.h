#ifndef  included_AMP_SubsetVectorUsingSuperDOFs
#define  included_AMP_SubsetVectorUsingSuperDOFs

#include "SubsetVector.h"
#include "VectorIndexer.h"

namespace AMP {
namespace LinearAlgebra {

  class SubsetVectorUsingSuperDOFs : public SubsetVector
  {
    protected:
      VectorIndexer::shared_ptr    d_Indexer;
      size_t         d_GlobalSize;
      size_t         d_LocalSize;
      size_t         d_LocalStartID;

    public:
      static  Vector::shared_ptr   view ( Vector::shared_ptr , Variable::shared_ptr );

      virtual void     setValuesByLocalID ( int , size_t * , const double * );
      virtual void     setLocalValuesByGlobalID ( int , size_t * , const double * );
      virtual void     addValuesByLocalID ( int , size_t * , const double * );
      virtual void     addLocalValuesByGlobalID ( int , size_t * , const double * );
      virtual void     getLocalValuesByGlobalID ( int , size_t * , double * ) const ;

      /* A SubsetVectorUsingSuperDOFs is assumed to be a subset of a vector on a volume and
         can interact with said volume.  For this reason, these methods return the size
         of the volume vector not the subet. */
      virtual size_t   getGlobalMaxID () const { return d_GlobalSize; }
      virtual size_t   getLocalMaxID () const { return d_LocalSize; }
      virtual size_t   getLocalStartID () const { return d_LocalStartID; }
  };

}
}

#endif
