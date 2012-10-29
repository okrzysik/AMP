
#ifndef included_AMP_MapOperatorParameters
#define included_AMP_MapOperatorParameters

#include "operators/OperatorParameters.h"
#include "ampmesh/Mesh.h"

namespace AMP {
namespace Operator {

class MapOperatorParameters : public OperatorParameters {
public :

      MapOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : OperatorParameters(db), d_MapComm(AMP_COMM_NULL) {  }

      virtual ~MapOperatorParameters() { }

      AMP_MPI d_MapComm;
      AMP::Mesh::Mesh::shared_ptr d_MapMesh;

};

}
}

#endif

