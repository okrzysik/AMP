
#ifndef included_AMP_MeshBasedOperatorParameters
#define included_AMP_MeshBasedOperatorParameters

#include "operators/OperatorParameters.h"
#include "ampmesh/Mesh.h"

namespace AMP {
  namespace Operator {
    class MeshBasedOperatorParameters : public OperatorParameters {
      public:
        MeshBasedOperatorParameters(const AMP::shared_ptr<AMP::Database> & db)
          : OperatorParameters(db) {  }

        virtual ~MeshBasedOperatorParameters() { }

        AMP::Mesh::Mesh::shared_ptr d_Mesh;
    };
  }
}


#endif


