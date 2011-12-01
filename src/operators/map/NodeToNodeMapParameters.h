#ifndef included_AMP_NodeToNodeMapParameters
#define included_AMP_NodeToNodeMapParameters

#include "operators/map/AsyncMapOperatorParameters.h"
#include "ampmesh/Mesh.h"

namespace AMP {
namespace Operator {

  class NodeToNodeMapParameters : public AMP::Operator::AsyncMapOperatorParameters
  {
    public:
      std::vector<AMP::Mesh::MeshElementID>         d_ids;
      std::vector<double>                           d_disps;
      std::map<size_t,size_t>                       d_RemoteToLocalId;
      size_t                                        d_NumPartners;

      NodeToNodeMapParameters ( const boost::shared_ptr<AMP::Database> &db )
        : AsyncMapOperatorParameters ( db )
      {
      }



  };


}
}

#endif
