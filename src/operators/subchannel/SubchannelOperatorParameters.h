
#ifndef included_AMP_SubchannelParameters
#define included_AMP_SubchannelParameters

#include "discretization/DOF_Manager.h"
#include "operators/OperatorParameters.h"
#include "operators/boundary/libmesh/RobinPhysicsModel.h"
#include "operators/subchannel/SubchannelPhysicsModel.h"

namespace AMP {
namespace Operator {

//! Parameter class to provide parameters to all subchannel classes
class SubchannelOperatorParameters : public OperatorParameters
{
public:
    //! Constructor
    explicit SubchannelOperatorParameters( const AMP::shared_ptr<AMP::Database> &db )
        : OperatorParameters( db ), d_initialize( false )
    {
    }

    //! Destructor
    virtual ~SubchannelOperatorParameters() {}

    // pointer to subchannel physics model
    AMP::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;

    AMP::shared_ptr<RobinPhysicsModel> d_dittusBoelterCoefficient;

    AMP::LinearAlgebra::Vector::shared_ptr d_frozenSolution;

    bool d_initialize; // Do we want to initialize the matrix

    // Clad properties
    std::vector<double> clad_x, clad_y, clad_d;
};
}
}

#endif
