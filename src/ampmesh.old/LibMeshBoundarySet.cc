#include "LibMeshBoundarySet.h"
#include "LibMeshAdapter.h"

#include "boundary_info.h"



namespace AMP { 
namespace Mesh {

  /*
  bool LibMeshBoundarySet::isMyFacet ( LibMeshAdapter::Facet &f , LibMeshAdapter &m )
  {
    for ( size_t i = 0 ; i != f.numNodes() ; i++ )
    {
      if ( m.getNode ( f.getNodeID ( i ) ).getNode().procID() > libMesh::processor_id() )
      {
        return false;
      }
    }
    return true;
  }
  */

  void LibMeshBoundarySet::build ( LibMeshAdapter &mesh )
  {
    std::vector<unsigned int>     n_ids;
    std::vector<unsigned int>     e_ids;
    std::vector<unsigned short int>     s_ids;
    std::vector<short int>        b_ids;
    mesh.getMesh().boundary_info->build_node_list_from_side_list ();
    mesh.getMesh().boundary_info->build_node_list ( n_ids , b_ids );
    AMP_MPI meshComm(AMP_COMM_WORLD);
    int rank = meshComm.getRank();
    int size = meshComm.getSize();
    unsigned int curFid = rank + 1u;
    for ( size_t i = 0 ; i != n_ids.size() ; i++ )
    {
      d_Boundaries[b_ids[i]].push_back ( &(mesh.getNode ( n_ids[i] ).getNode() ) );
    }

    b_ids.clear();
    mesh.getMesh().boundary_info->build_side_list ( e_ids , s_ids , b_ids );
    for ( size_t i = 0 ; i != e_ids.size() ; i++ )
    {
      LibMeshAdapter::Element curElem = mesh.getElement ( e_ids[i] );
      LibMeshAdapter::Facet curFacet = curElem.getSide ( s_ids[i] );
//   There are two algorithms for doing this.  The first is to check the element the boundary
//   is on.  The second is to check the nodes that comprise the facet.
      if ( curElem.isOwned() )
//    if ( isMyFacet ( curFacet , mesh ) )
      {
        ::Elem *cur_side = curFacet.getPtr().release();
        cur_side->set_id ( curFid );
        d_ElementsOfFacets[b_ids[i]][curFid] = &(curElem.getElem());
        curFid += size;
        d_SideBoundaries[b_ids[i]].push_back ( cur_side );
        d_SidesToDelete.insert ( cur_side );
      }
    }
  }

  LibMeshBoundarySet::~LibMeshBoundarySet ()
  {
    std::set< ::Elem * >::iterator cur_side = d_SidesToDelete.begin();
    while ( cur_side != d_SidesToDelete.end() )
    {
      delete *cur_side;
      cur_side++;
    }
  }


  LibMeshBoundarySet::NodeListIterator  LibMeshBoundarySet::beginBoundary ( short int which ) 
  {
    return d_Boundaries[which].begin();
  }

  LibMeshBoundarySet::NodeListIterator  LibMeshBoundarySet::endBoundary ( short int which ) 
  {
    return d_Boundaries[which].end();
  }

  LibMeshBoundarySet::SideListIterator  LibMeshBoundarySet::beginSideBoundary ( short int which ) 
  {
    return d_SideBoundaries[which].begin();
  }

  LibMeshBoundarySet::SideListIterator  LibMeshBoundarySet::endSideBoundary ( short int which ) 
  {
    return d_SideBoundaries[which].end();
  }

}
}

