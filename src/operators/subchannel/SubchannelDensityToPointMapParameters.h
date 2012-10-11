#ifndef included_AMP_SubchannelDensityToPointMapParameters
#define included_AMP_SubchannelDensityToPointMapParameters

#include "operators/OperatorParameters.h"
#include "utils/AMP_MPI.h"
#include "operators/subchannel/SubchannelPhysicsModel.h"


namespace AMP {
namespace Operator {


/**
 * \class SubchannelDensityToPointMapParameters
 * \brief Parameters for the SubchannelDensityToPointMap
 */
class SubchannelDensityToPointMapParameters : public AMP::Operator::OperatorParameters
{
public :

    //! Default constructors
    SubchannelDensityToPointMapParameters(): 
        OperatorParameters(boost::shared_ptr<AMP::Database>()), 
        d_comm(AMP_COMM_WORLD) {};

    //! Deconstructor
    virtual ~SubchannelDensityToPointMapParameters() { }

    //! Comm over which the points will be gathered (default is AMP_COMM_WORLD)
    AMP_MPI d_comm;

    // List of the coordinates of the points
    std::vector<double> x;      //!< x-coordinate of the points to fill
    std::vector<double> y;      //!< y-coordinate of the points to fill
    std::vector<double> z;      //!< z-coordinate of the points to fill
    
    // Subchannel physics model
    boost::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;
};


}
}

#endif

