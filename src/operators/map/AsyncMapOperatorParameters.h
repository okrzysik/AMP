#ifndef included_AMP_AsyncMapOperatorParameters
#define included_AMP_AsyncMapOperatorParameters

#include "operators/AsynchronousOperatorParameters.h"
#include "utils/AMP_MPI.h"

#include "ampmesh/MeshManager.h"

namespace AMP {
namespace Operator {

  class AsyncMapOperatorParameters : public AsynchronousOperatorParameters
  {
    public:
      typedef  AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  BNIterator;

      AMP_MPI          d_MapComm;
      BNIterator       d_SetBegin;
      BNIterator       d_SetEnd;

      bool             d_IsMaster;
      int              d_ToMasterCommTag;
      int              d_ToSlaveCommTag;

      AsyncMapOperatorParameters ( const boost::shared_ptr<AMP::Database> &db );

      virtual ~AsyncMapOperatorParameters ();

  };

}
}

#endif

