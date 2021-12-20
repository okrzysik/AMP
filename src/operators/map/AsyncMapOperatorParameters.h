#ifndef included_AMP_AsyncMapOperatorParameters
#define included_AMP_AsyncMapOperatorParameters

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/operators/AsynchronousOperatorParameters.h"
#include "AMP/utils/AMP_MPI.h"


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

    explicit AsyncMapOperatorParameters( std::shared_ptr<AMP::Database> db );

    virtual ~AsyncMapOperatorParameters();
};
} // namespace Operator
} // namespace AMP

#endif
