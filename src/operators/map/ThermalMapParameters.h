#ifndef included_AMP_ThermalMapParameters
#define included_AMP_ThermalMapParameters

#include "operators/map/AsyncMapOperatorParameters.h"
#include "ampmesh/Mesh.h"

namespace AMP {
namespace Operator {

  class ThermalMapParameters : public AsyncMapOperatorParameters
  {
    public:
      std::vector<int>          d_ids;
      std::vector<double>       d_disps;
      std::map<size_t,size_t>   d_RemoteToLocalId;
      size_t                    d_NumPartners;

      ThermalMapParameters ( const boost::shared_ptr<AMP::Database> &db )
        : AsyncMapOperatorParameters ( db )
      {
      }



  };


}
}

#endif
