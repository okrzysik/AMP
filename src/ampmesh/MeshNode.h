#ifndef included_AMP_MeshNode
#define included_AMP_MeshNode

#include "MeshObject.h"

namespace AMP { 
namespace Mesh {

  class MeshNode : public MeshObject
                 , public MeshPoint
  {
    protected:
      MeshNode () {}

    public:
      virtual ~MeshNode () {}

      virtual void translate ( double x , double y , double z ) = 0;
  };


}
}
#endif
