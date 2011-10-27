#ifndef included_AMP_MeshCollectionVariable_h
#define included_AMP_MeshCollectionVariable_h

#include "MeshVariable.h"

namespace AMP { 
namespace Mesh {

  // These are stubs for use when we incorporate more operators.


  class MeshCollectionVariable : public AMP::LinearAlgebra::Variable
  {
    private:
      AMP::LinearAlgebra::Variable::shared_ptr    d_BaseVariable;

    public:
      MeshCollectionVariable ( AMP::LinearAlgebra::Variable::shared_ptr base ) :  Variable ( base->getName() ) ,
                                                              d_BaseVariable ( base ) {}

      virtual bool operator == ( const Variable &rhs ) const
      {
        if ( rhs.isA<MeshVariable<MESH_TYPE> > () )
          return *d_BaseVariable == rhs;
        if ( rhs.isA<MeshCollectionVariable> () )
          return *d_BaseVariable == *(rhs.castTo<MeshCollectionVariable>().d_BaseVariable);
        return false;
      }

  };

}
}

#endif
