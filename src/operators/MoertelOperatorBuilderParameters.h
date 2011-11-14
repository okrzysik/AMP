#ifndef included_MoertelOperatorBuilderParameters_h
#define included_MoertelOperatorBuilderParameters_h

#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "utils/AMP_MPI.h"

namespace AMP {
namespace Operator {

  class MoertelOperatorBuilderParameters : public ParameterBase
  {
    public:
      AMP::Mesh::Mesh::shared_ptr                   d_MasterVolume;
      AMP::Mesh::Mesh::shared_ptr                       d_SlaveVolume;
      AMP::Discretization::DOFManager::shared_ptr       d_MasterDOFMap;
      AMP::Discretization::DOFManager::shared_ptr       d_SlaveDOFMap;
      short int                                         d_MasterSurface;
      short int                                         d_SlaveSurface;
      AMP_MPI                                           d_Comm;
      boost::shared_ptr<Database>                       d_DB;
  };



}
}

#endif
