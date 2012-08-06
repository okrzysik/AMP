
#ifndef included_AMP_StructuredMeshHelper
#define included_AMP_StructuredMeshHelper

#include "utils/Utilities.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElementVectorIterator.h"

namespace AMP{
namespace Mesh{

  class StructuredMeshHelper 
  {
    public : 

      static AMP::Mesh::MeshIterator getXYFaceIterator(AMP::Mesh::Mesh::shared_ptr subChannel, int ghostWidth);

  };

}
}
#endif
