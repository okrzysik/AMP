#include "MeshCollectionVector.h"

namespace AMP { 
namespace Mesh {

  AMP::LinearAlgebra::Vector::shared_ptr MeshCollectionVector::cloneVector(const AMP::LinearAlgebra::Variable::shared_ptr name) const
  {
    AMP::LinearAlgebra::Vector::shared_ptr  retVal = AMP::LinearAlgebra::Vector::shared_ptr ( new MeshCollectionVector ( name ) );
    MeshCollectionVector  &ret = retVal->castTo<MeshCollectionVector> ();
    ret.d_vVectors.resize ( d_vVectors.size() );
    ret.d_vGlobalOffsets.resize ( d_vGlobalOffsets.size() );
    ret.d_vLocalOffsets.resize ( d_vLocalOffsets.size() );
    for ( size_t i = 0 ; i != d_vGlobalOffsets.size() ; i++ )
    {
      ret.d_vGlobalOffsets[i] = d_vGlobalOffsets[i];
      ret.d_vLocalOffsets[i] = d_vLocalOffsets[i];
    }
    for ( size_t i = 0 ; i != d_vVectors.size() ; i++ )
    {
      ret.d_vVectors[i] = d_vVectors[i]->cloneVector ();
    }
    return retVal;
  }

}
}

