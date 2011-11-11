#ifndef included_ElementOperationFactory
#define  included_ElementOperationFactory

/* Boost files */
#include "boost/shared_ptr.hpp"

#include "operators/ElementOperation.h"


namespace AMP {
namespace Operator {

class ElementOperationFactory{
 public:
  ElementOperationFactory(){}
  ~ElementOperationFactory(){}

  static boost::shared_ptr<ElementOperation> createElementOperation(boost::shared_ptr<AMP::Database>  input_db);

};  

}
}

#endif
 
