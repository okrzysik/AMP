#ifndef included_ParameterFactory
#define included_ParameterFactory

/* AMP files */
#include "ampmesh/MeshManager.h"
#include "utils/Database.h"
#include "operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

class ParameterFactory{

 public:
  ParameterFactory(){}
  ~ParameterFactory(){}

  static boost::shared_ptr<OperatorParameters>  createParameter(boost::shared_ptr<AMP::Database>  input_db, AMP::Mesh::MeshManager::Adapter::shared_ptr  mesh);

};
  
}
}

#endif
