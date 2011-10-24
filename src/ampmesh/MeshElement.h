#ifndef included_AMP_MeshElement
#define included_AMP_MeshElement

#include "MeshObject.h"

namespace AMP { 
namespace Mesh {

  class MeshElement : public MeshObject
  {
    protected:
      MeshElement () {}

    public:
      virtual ~MeshElement () {}

      virtual size_t   numNodes () const = 0;
      virtual size_t   numSides () const = 0;

      virtual bool     hasNeighbor ( size_t k ) const = 0;
      virtual size_t   getNodeID ( size_t k ) const = 0;

      virtual bool    isOwned () const = 0;

      virtual double  volume () const = 0;
  };

}
}
#endif
