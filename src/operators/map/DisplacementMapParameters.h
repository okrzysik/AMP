#ifndef included_AMP_DisplacementMapParameters
#define included_AMP_DisplacementMapParameters

#include "operators/map/AsyncMapOperatorParameters.h"
#include "ampmesh/MeshManager.h"

namespace AMP {
namespace Operator {

  class DisplacementMapParameters : public AMP::Operator::AsyncMapOperatorParameters
  {
    public:
      std::vector<int>          d_ids;
      std::vector<double>       d_disps;
      std::map<size_t,size_t>   d_RemoteToLocalId;
      size_t                    d_NumPartners;

      DisplacementMapParameters ( const boost::shared_ptr<AMP::Database> &db )
        : AMP::Operator::AsyncMapOperatorParameters ( db )
      {
      }



  };


}
}

#endif
