#ifndef included_AMP_LibMeshNodeElemMap_h
#define included_AMP_LibMeshNodeElemMap_h

#include <map>
#include <vector>
#include <cstddef>

class Elem;

namespace AMP { 
namespace Mesh {

  class LibMeshAdapter;

  class LibMeshNodeElemMap
  {
    public:
      typedef   std::vector< ::Elem *>               ElemList;
      typedef   std::map<size_t , ElemList>          NodeList;
      typedef   NodeList::iterator                   NodeListIterator;
      typedef   ElemList::iterator                   ElemListIterator;

    private:
      NodeList     d_Nodes;

    public:
      LibMeshNodeElemMap () {}

      void              build ( LibMeshAdapter & );

      ElemListIterator  beginElement ( size_t which );
      ElemListIterator  endElement   ( size_t which );
  };

}
}

#endif
