#include "AsyncMapOperator.h"
#include "AsyncMapOperatorParameters.h"

namespace AMP {
namespace Operator {

size_t mapTagOffset = 0;

  AsyncMapOperator::AsyncMapOperator ( const boost::shared_ptr <OperatorParameters> &p )
    : AsynchronousOperator ( p )
  {
    boost::shared_ptr<AsyncMapOperatorParameters> params = boost::dynamic_pointer_cast<AsyncMapOperatorParameters>(p);
    d_IsMaster = params->d_IsMaster;
  }

  AsyncMapOperator::~AsyncMapOperator ()
  {
  }

}
}

