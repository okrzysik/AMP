#ifndef included_AMP_LibMeshPoint_h
#define included_AMP_LibMeshPoint_h

#include "MeshPoint.h"

#include "mesh.h"
#include "elem.h"
#include "node.h"

namespace AMP { 
namespace Mesh {

  class LibMeshPoint : public MeshPoint
  {
    private:
      ::Point   &d_Pt;

    public:
      LibMeshPoint ( ::Point &rhs );
      virtual ~LibMeshPoint ();

      virtual  double    x () const;
      virtual  double    y () const;
      virtual  double    z () const;

      virtual  double   &x ();
      virtual  double   &y ();
      virtual  double   &z ();

      virtual  double  &operator ()( size_t i );
      virtual  double   operator ()( size_t i ) const;

  };

}
}

#include "LibMeshPoint.inline.h"
#endif
