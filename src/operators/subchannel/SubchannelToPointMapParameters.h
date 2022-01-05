#ifndef included_AMP_SubchannelToPointMapParameters
#define included_AMP_SubchannelToPointMapParameters

#include "AMP/operators/OperatorParameters.h"
#include "AMP/operators/subchannel/SubchannelPhysicsModel.h"
#include "AMP/utils/AMP_MPI.h"


namespace AMP::Operator {


/**
 * \class SubchannelToPointMapParameters
 * \brief Parameters for the SubchannelToPointMap
 */
class SubchannelToPointMapParameters : public AMP::Operator::OperatorParameters
{
public:
    //! Default constructors
    SubchannelToPointMapParameters()
        : OperatorParameters( std::shared_ptr<AMP::Database>() ), d_comm( AMP_COMM_WORLD ){};

    //! Deconstructor
    virtual ~SubchannelToPointMapParameters() {}

    //! Comm over which the points will be gathered (default is AMP_COMM_WORLD)
    AMP_MPI d_comm;

    // List of the coordinates of the points
    std::vector<double> x; //!< x-coordinate of the points to fill
    std::vector<double> y; //!< y-coordinate of the points to fill
    std::vector<double> z; //!< z-coordinate of the points to fill

    // Subchannel physics model
    std::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;

    // Output variable (may be null on processors where x is empty)
    // Valid variables are: "Density", "Temperature"
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_outputVar;
};
} // namespace AMP::Operator

#endif
