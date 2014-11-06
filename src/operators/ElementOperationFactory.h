#ifndef included_ElementOperationFactory
#define  included_ElementOperationFactory

/* Boost files */
#include "utils/shared_ptr.h"

#include "operators/ElementOperation.h"


namespace AMP {
namespace Operator {

class ElementOperationFactory{
 public:
  ElementOperationFactory(){}
  ~ElementOperationFactory(){}

  static AMP::shared_ptr<ElementOperation> createElementOperation(AMP::shared_ptr<AMP::Database>  input_db);

};  

}
}

#endif
 
