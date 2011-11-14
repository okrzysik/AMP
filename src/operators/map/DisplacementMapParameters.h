#ifndef included_AMP_DisplacementMapParameters
#define included_AMP_DisplacementMapParameters

#include "operators/map/AsyncMapOperatorParameters.h"

namespace AMP {
namespace Operator {

  class DisplacementMapParameters : public AMP::Operator::AsyncMapOperatorParameters
  {
    public:
      std::vector<Mesh::MeshElementID>          d_ids;              //!< Ids of the mesh nodes
      std::vector<double>                       d_disps;
      std::map<size_t,size_t>                   d_RemoteToLocalId;
      size_t                                    d_NumPartners;

      DisplacementMapParameters ( const boost::shared_ptr<AMP::Database> &db )
        : AMP::Operator::AsyncMapOperatorParameters ( db )
      {
      }



  };


}
}

#endif
