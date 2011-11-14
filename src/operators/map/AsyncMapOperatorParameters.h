#ifndef included_AMP_AsyncMapOperatorParameters
#define included_AMP_AsyncMapOperatorParameters

#include "operators/AsynchronousOperatorParameters.h"
#include "utils/AMP_MPI.h"


namespace AMP {
namespace Operator {

  class AsyncMapOperatorParameters : public AsynchronousOperatorParameters
  {
    public:

      AMP_MPI                       d_MapComm;
      AMP::Mesh::MeshIterator       d_BoundaryNodeIterator;

      bool             d_IsMaster;
      int              d_ToMasterCommTag;
      int              d_ToSlaveCommTag;

      AsyncMapOperatorParameters ( const boost::shared_ptr<AMP::Database> &db );

      virtual ~AsyncMapOperatorParameters ();

  };

}
}

#endif

