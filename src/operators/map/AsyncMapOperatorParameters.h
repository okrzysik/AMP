#ifndef included_AMP_AsyncMapOperatorParameters
#define included_AMP_AsyncMapOperatorParameters

#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "operators/AsynchronousOperatorParameters.h"
#include "utils/AMP_MPI.h"


namespace AMP {
namespace Operator {

class AsyncMapOperatorParameters : public AsynchronousOperatorParameters
{
public:
    AMP_MPI d_MapComm;
    AMP::Mesh::Mesh::shared_ptr d_Mesh1;
    AMP::Mesh::Mesh::shared_ptr d_Mesh2;
    int d_BoundaryID1;
    int d_BoundaryID2;
    int d_commTag;
    bool callMakeConsistentSet;

    explicit AsyncMapOperatorParameters( const AMP::shared_ptr<AMP::Database> &db );

    virtual ~AsyncMapOperatorParameters();
};
}
}

#endif
