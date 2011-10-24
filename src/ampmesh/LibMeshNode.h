#ifndef included_AMP_LibMeshNode_h
#define included_AMP_LibMeshNode_h

#include "MeshNode.h"
#include "MeshElement.h"

#include "mesh.h"
#include "elem.h"
#include "node.h"

namespace AMP { 
namespace Mesh {

  class LibMeshNode : public MeshNode
  {
    private:
      ::Node  *d_Node;

    public:
      // This needs to be fixed at some point.  No null ptrs!
      LibMeshNode ();

      LibMeshNode ( ::Node *n );
      LibMeshNode ( const LibMeshNode &rhs );

            ::Node &getNode();
      const ::Node &getNode() const;

      size_t   globalID () const;

      virtual void translate ( double x , double y , double z );

      virtual double   &x ();
      virtual double   &y ();
      virtual double   &z ();

      virtual double    x () const;
      virtual double    y () const;
      virtual double    z () const;


      virtual double  operator () ( size_t t ) const;
      virtual double &operator () ( size_t t );

      virtual bool    isOwned () const;
      virtual size_t  procID () const;
  };


}
}

#include "LibMeshNode.inline.h"
#endif
