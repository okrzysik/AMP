#ifndef included_AMP_LibMeshBoundarySet_h
#define included_AMP_LibMeshBoundarySet_h

#include <map>
#include <vector>



// Boundary ID 1 = top
// Boundary ID 2 = bottom
// Boundary ID 4 = outside
// Boundary ID 8 = inside

// Boundary ID 2048 = Dirichlet face  (for when you don't have contact to make problems well-posed)
// Boundary ID 2049 = Dirichlet face top
// Boundary ID 2050 = Dirichlet face bottom
// Boundary ID 2052 = Dirichlet face outside
// Boundary ID 2056 = Dirichlet face inside

// This is the LibMesh node
#include "elem.h"
#include "node.h"

namespace AMP { 
namespace Mesh {

  class LibMeshAdapter;

  class LibMeshBoundarySet 
  {
    public:
      typedef   std::vector< ::Node *>               NodeList;
      typedef   std::vector< ::Elem * >              ElemList;
      typedef   std::map<short int,NodeList>         BoundaryList;
      typedef   std::map<short int,ElemList>         SideBoundaryList;
      typedef   BoundaryList::iterator               BoundaryListIterator;
      typedef   NodeList::iterator                   NodeListIterator;
      typedef   ElemList::iterator                   SideListIterator;

    private:
      BoundaryList                   d_Boundaries;
      SideBoundaryList               d_SideBoundaries;
      std::set< ::Elem * >           d_SidesToDelete;
      std::map<short int , std::map<size_t , ::Elem *> >   d_ElementsOfFacets;

  //    bool  isMyFacet ( LibMeshAdapter::Facet &f , LibMeshAdapter &m );

    public:
      LibMeshBoundarySet () {}
     ~LibMeshBoundarySet ();

      void              build ( LibMeshAdapter & );

      NodeListIterator  beginBoundary ( short int which );
      NodeListIterator  endBoundary ( short int which );

      SideListIterator  beginSideBoundary ( short int which );
      SideListIterator  endSideBoundary ( short int which );

      ::Elem *          getElementOfFacet ( short int bid , size_t fid )
      {
        return d_ElementsOfFacets[bid][fid];
      }
  };

}
}

#endif
