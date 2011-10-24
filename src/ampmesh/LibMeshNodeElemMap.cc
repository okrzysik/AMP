#include "LibMeshNodeElemMap.h"
#include "LibMeshAdapter.h"


namespace AMP { 
namespace Mesh {


void LibMeshNodeElemMap::build ( LibMeshAdapter &mesh )
{
    LibMeshAdapter::ElementIterator  cur_elem = mesh.beginElement();
    LibMeshAdapter::ElementIterator  end_elem = mesh.endElement();
    while ( cur_elem != end_elem ) {
        ::Elem *elem_ptr = &(cur_elem->getElem());
        for ( size_t i = 0 ; i != cur_elem->numNodes() ; i++ )
            d_Nodes[cur_elem->getNodeID(i)].push_back( elem_ptr );
        cur_elem++;
    }
}


LibMeshNodeElemMap::ElemListIterator  LibMeshNodeElemMap::beginElement ( size_t which )
{
    return d_Nodes[which].begin();
}


LibMeshNodeElemMap::ElemListIterator  LibMeshNodeElemMap::endElement ( size_t which )
{
    return d_Nodes[which].end();
}


}
}
