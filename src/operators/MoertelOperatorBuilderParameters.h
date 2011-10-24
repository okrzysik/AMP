#ifndef included_MoertelOperatorBuilderParameters_h
#define included_MoertelOperatorBuilderParameters_h

#include "ampmesh/MeshManager.h"
#include "utils/AMP_MPI.h"

namespace AMP {
namespace Operator {

  class MoertelOperatorBuilderParameters : public ParameterBase
  {
    public:
      AMP::Mesh::MeshManager::Adapter::shared_ptr   d_MasterVolume;
      AMP::Mesh::MeshManager::Adapter::shared_ptr   d_SlaveVolume;
      AMP::Mesh::DOFMap::shared_ptr                 d_MasterDOFMap;
      AMP::Mesh::DOFMap::shared_ptr                 d_SlaveDOFMap;
      short int                          d_MasterSurface;
      short int                          d_SlaveSurface;
      AMP_MPI                            d_Comm;
      boost::shared_ptr<Database>        d_DB;
  };



}
}

#endif
