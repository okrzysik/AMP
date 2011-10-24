
#ifndef included_AMP_MapOperatorParameters
#define included_AMP_MapOperatorParameters

#include "operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

  class MapOperatorParameters : public OperatorParameters {
    public :

      MapOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : OperatorParameters(db) {  }

      ~MapOperatorParameters() { }

      AMP::Mesh::MeshManager::Adapter::shared_ptr d_MapAdapter;

    protected :

    private :

  };

}
}

#endif

