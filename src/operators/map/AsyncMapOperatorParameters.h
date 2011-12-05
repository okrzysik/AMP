#ifndef included_AMP_AsyncMapOperatorParameters
#define included_AMP_AsyncMapOperatorParameters

#include "operators/AsynchronousOperatorParameters.h"
#include "ampmesh/Mesh.h"
#include "utils/AMP_MPI.h"
#include "discretization/DOF_Manager.h"


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
      AMP::Discretization::DOFManager::shared_ptr  d_DOFManager;
      int                           d_commTag;

      AsyncMapOperatorParameters ( const boost::shared_ptr<AMP::Database> &db );

      virtual ~AsyncMapOperatorParameters ();

  };

}
}

#endif

