#ifndef included_AMP_LibMeshElement
#define included_AMP_LibMeshElement

#include "MeshNode.h"
#include "MeshElement.h"

#include "mesh.h"
#include "elem.h"
#include "node.h"

namespace AMP { 
namespace Mesh {

  typedef ::ElemType ElementType;

  template <typename PTR_TYPE>
  class LibMeshElementTmpl : public MeshElement
  {
    private:
      PTR_TYPE   d_Elem;
      ::Elem    *d_Parent;

    public:
      // This needs to be fixed at some point.  No null ptrs!
      LibMeshElementTmpl ();

      LibMeshElementTmpl ( PTR_TYPE e );

      // This breaks auto_ptr!
      LibMeshElementTmpl ( const LibMeshElementTmpl &rhs ) ;

      PTR_TYPE  getPtr ();

      // Gotta love shared pointers.
      void nonConstConstruct ( const LibMeshElementTmpl rhs );


      void  setParent ( ::Elem *e ) { d_Parent = e; }
      ::Elem *getParent () { return d_Parent; }
      virtual size_t   globalID () const;
      virtual bool  operator < ( const MeshObject &rhs ) const;
      virtual bool  operator == ( const LibMeshElementTmpl &rhs ) const;

            ::Elem &getElem ();
      const ::Elem &getElem () const;

      virtual size_t   numNodes () const;
      virtual size_t   numSides () const;

      virtual bool     hasNeighbor ( size_t k ) const;
      LibMeshElementTmpl< AutoPtr < ::Elem> >  getSide ( size_t k , bool proxy = false );

      LibMeshPoint  getPoint ( size_t k );
      virtual size_t   getNodeID ( size_t k ) const;

      ElementType   getType () const;
      virtual bool    isOwned () const;
      size_t          procID () const;

      double volume() const;
  };

  typedef LibMeshElementTmpl< ::Elem *>            LibMeshElement;
  typedef LibMeshElementTmpl< AutoPtr< ::Elem> >   LibMeshSide;


}
}

#include "LibMeshElement.tmpl.h"
#endif
