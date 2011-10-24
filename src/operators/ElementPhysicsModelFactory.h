#ifndef included_AMP_ElementPhysicsModelFactory
#define included_AMP_ElementPhysicsModelFactory

/* Boost files */
#include "boost/shared_ptr.hpp"

#include "ElementPhysicsModel.h"
#include "operators/boundary/RobinPhysicsModel.h"

namespace AMP {
namespace Operator {

class ElementPhysicsModelFactory{
 public:
  ElementPhysicsModelFactory(){}
  ~ElementPhysicsModelFactory(){}

  static boost::shared_ptr<ElementPhysicsModel> createElementPhysicsModel(boost::shared_ptr<AMP::Database>  input_db);

};
  
}
}

#endif
