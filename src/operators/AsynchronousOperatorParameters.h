#ifndef included_AMP_AsynchronousOperatorParameters
#define included_AMP_AsynchronousOperatorParameters

#include "operators/OperatorParameters.h"

namespace AMP {
namespace Operator {


  class AsynchronousOperatorParameters : public OperatorParameters
  {
    public:
      char    d_AsynchronousConstructionParam;
      size_t  d_ConstructionPhase;

      AsynchronousOperatorParameters ( const boost::shared_ptr<AMP::Database> &db )
        : OperatorParameters ( db )
      {
      }



  };


}
}

#endif
