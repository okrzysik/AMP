#ifndef included_AMP_AsyncMapOperatorParameters
#define included_AMP_AsyncMapOperatorParameters

#include "operators/AsynchronousOperatorParameters.h"
#include "ampmesh/Mesh.h"
#include "utils/AMP_MPI.h"


namespace AMP {
namespace Operator {

  class AsyncMapOperatorParameters : public AsynchronousOperatorParameters
  {
    public:

      AMP_MPI                       d_MapComm;
      AMP::Mesh::Mesh::shared_ptr   d_Mesh1;
      AMP::Mesh::Mesh::shared_ptr   d_Mesh2;
      int                           d_BoundaryID1;
      int                           d_BoundaryID2;

      bool             d_IsMaster;
      int              d_ToMasterCommTag;
      int              d_ToSlaveCommTag;

      AsyncMapOperatorParameters ( const boost::shared_ptr<AMP::Database> &db );

      virtual ~AsyncMapOperatorParameters ();

  };

}
}

#endif

