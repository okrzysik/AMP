#ifndef included_AMP_MeshObject_h
#define included_AMP_MeshObject_h

#include <boost/shared_ptr.hpp>
#include "utils/Castable.h"

namespace AMP { 
namespace Mesh {

  class MeshObject : public Castable
  {
    public:
      typedef boost::shared_ptr<MeshObject>    shared_ptr;
      MeshObject () {}
      virtual ~MeshObject () {}

      virtual size_t      globalID() const = 0;
      virtual bool        isOwned () const = 0;
      virtual size_t      procID () const = 0;

      virtual bool        operator < ( const MeshObject &obj ) const
      {
        return globalID() < obj.globalID();
      }
  };

}
}

#endif
