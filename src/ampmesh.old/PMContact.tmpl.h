
#include "MeshUtils.h"
#include "math.h"

namespace AMP { 
namespace Mesh {

  template <typename MANAGER>
  double  PMContact<MANAGER>::computeShiftUp ( MeshManagerPtr  manager , std::vector<std::string> &prefix)
  {
    std::string  mesh_name = manager->beginMeshNames()->first;
    int color = 0;
    for ( size_t i = 0 ; i != prefix.size() ; i++ )
    {
      std::string  mesh_prefix = mesh_name.substr ( 0 , prefix[i].size() );
      color = mesh_prefix == prefix[i] ? 1 : 0;
    }
    AMP_MPI globalComm(AMP_COMM_WORLD);
    int rank = globalComm.getRank();
    AMP_MPI semiLocalComm = globalComm.split(color,rank);
    if ( color )
    {
      MeshAdapterPtr  meshAdapter = *(manager->beginMeshes());
      min_max_struct<simple_point> boundingBox = computeExtremeCoordinates<MeshAdapter> ( meshAdapter );
      double meshHeight = boundingBox.max.z - boundingBox.min.z;
      int meshRank = (manager->getMeshComm()).getRank();
      if ( meshRank != 0 )
      {
        meshHeight = 0.0;
      }
      double totalHeight;
      semiLocalComm.sumScan(&meshHeight,&totalHeight,1);
      double bottomHeight = totalHeight - meshHeight;
      return bottomHeight;
    }
    return 0.0;
  }

  template <typename MANAGER>
  void  PMContact<MANAGER>::restrictOuterRadius ( MeshManagerPtr manager , double r , const std::string &prefix )
  {
    std::string  mesh_name = manager->beginMeshNames()->first;
    std::string  mesh_prefix = mesh_name.substr ( 0 , prefix.size() );
    if ( mesh_prefix == prefix )
    {
      MeshAdapterPtr  meshAdapter = *(manager->beginMeshes());
      typename MeshAdapter::NodeIterator  curNode = meshAdapter->beginUsedNode ();
      while ( curNode != meshAdapter->endUsedNode() )
      {
        double &cur_x = curNode->x();
        double &cur_y = curNode->y();
        double  cur_r = sqrt ( cur_x*cur_x + cur_y*cur_y );
        if ( cur_r > r )
        {
          double cur_theta = atan2 ( cur_y , cur_x );
          cur_x = r * cos ( cur_theta );
          cur_y = r * sin ( cur_theta );
        }
        curNode++;
      }
    }
  }
  template <typename MANAGER>
  void  PMContact<MANAGER>::restrictInnerRadius ( MeshManagerPtr manager , double r , const std::string &prefix )
  {
    std::string  mesh_name = manager->beginMeshNames()->first;
    std::string  mesh_prefix = mesh_name.substr ( 0 , prefix.size() );
    if ( mesh_prefix == prefix )
    {
      MeshAdapterPtr  meshAdapter = *(manager->beginMeshes());
      typename MeshAdapter::NodeIterator  curNode = meshAdapter->beginUsedNode ();
      while ( curNode != meshAdapter->endUsedNode() )
      {
        double &cur_x = curNode->x();
        double &cur_y = curNode->y();
        double  cur_r = sqrt ( cur_x*cur_x + cur_y*cur_y );
        if ( cur_r < r )
        {
          double cur_theta = atan2 ( cur_y,cur_x );
          cur_x = r * cos ( cur_theta );
          cur_y = r * sin ( cur_theta );
        }
        curNode++;
      }
    }
  }

}
}


