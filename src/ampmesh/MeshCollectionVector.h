#ifndef included_AMP_MeshCollectionVector
#define included_AMP_MeshCollecitonVector

#include "vectors/MultiVector.h"
#include "MeshCollectionVariable.h"
#include "utils/AMP_MPI.h"

namespace AMP { 
namespace Mesh {

  class MeshCollectionVector : public AMP::LinearAlgebra::MultiVector
  {
    protected:
      MeshCollectionVector ( AMP::LinearAlgebra::Variable::shared_ptr n ) : MultiVector ( n ) {}

    public:
      template <typename MESHTYPE>
      static AMP::LinearAlgebra::Vector::shared_ptr    create ( AMP::LinearAlgebra::Variable::shared_ptr name , AMP_MPI comm );

      virtual std::string type() const { return "Mesh Collection Vector"; }
      virtual AMP::LinearAlgebra::Vector::shared_ptr cloneVector(const AMP::LinearAlgebra::Variable::shared_ptr name) const;
      virtual void  addVector ( AMP::LinearAlgebra::Vector::shared_ptr v );
  };

  inline
  void MeshCollectionVector::addVector ( AMP::LinearAlgebra::Vector::shared_ptr v )
  {
//    AMP_ASSERT ( getVariable()->isSameTypeAs ( v->getVariable() ) );
//    AMP_ASSERT ( getVariable()->getName() == v->getVariable()->getName() );
    AMP::LinearAlgebra::MultiVector::addVector ( v );
  }

  template <typename MESHTYPE>
  AMP::LinearAlgebra::Vector::shared_ptr    MeshCollectionVector::create ( AMP::LinearAlgebra::Variable::shared_ptr name , AMP_MPI comm )
  {
    AMP::LinearAlgebra::Vector::shared_ptr  retval;
    if ( name->isA<MeshVariable<MESHTYPE> > () )
    {
      AMP::LinearAlgebra::Variable::shared_ptr  newVar ( new MeshCollectionVariable ( name ) );
      retval = AMP::LinearAlgebra::Vector::shared_ptr ( new MeshCollectionVector ( newVar ) );
    }
    else if ( name->isA<MeshCollectionVariable> () )
    {
      retval = AMP::LinearAlgebra::Vector::shared_ptr ( new MeshCollectionVector ( name ) );
    }
    if ( !retval )
      AMP_ERROR( "Invalid variable type for MeshCollectionVector::create" );
    retval->castTo<MeshCollectionVector>().d_Comm = comm.getCommunicator();
    return retval;
  }

}
}

#endif


