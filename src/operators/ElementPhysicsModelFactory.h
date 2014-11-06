#ifndef included_AMP_ElementPhysicsModelFactory
#define included_AMP_ElementPhysicsModelFactory

/* Boost files */
#include "utils/shared_ptr.h"

#include "operators/ElementPhysicsModel.h"

namespace AMP {
namespace Operator {

class ElementPhysicsModelFactory{
 public:
  ElementPhysicsModelFactory(){}
  ~ElementPhysicsModelFactory(){}

  static AMP::shared_ptr<ElementPhysicsModel> createElementPhysicsModel(AMP::shared_ptr<AMP::Database>  input_db);

};
  
}
}

#endif
