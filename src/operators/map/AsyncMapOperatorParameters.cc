

#include "AsyncMapOperatorParameters.h"

namespace AMP {
namespace Operator {

  AsyncMapOperatorParameters::AsyncMapOperatorParameters ( const boost::shared_ptr<AMP::Database> &db )
    : AsynchronousOperatorParameters ( db )
  {
  }

  AsyncMapOperatorParameters::~AsyncMapOperatorParameters ()
  {
  }

}
}
