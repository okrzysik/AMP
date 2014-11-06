#ifndef included_AMP_SubchannelToPointMapMap
#define included_AMP_SubchannelToPointMapMap

#include "utils/AMP_MPI.h"
#include "ampmesh/Mesh.h"
#include "operators/Operator.h"
#include "operators/subchannel/SubchannelToPointMapParameters.h"
#include "operators/subchannel/SubchannelPhysicsModel.h"

namespace AMP {
namespace Operator {


/**
 * \class SubchannelToPointMap
 * \brief A class used to map subchannel properties to points
 * \details  This class provides routines for mapping the subchannel flow  
 *     properties (Density, Temperature) to a set of points provided.
 */
class SubchannelToPointMap : public AMP::Operator::Operator
{
public :

    //! Default constructor
    SubchannelToPointMap(const AMP::shared_ptr<SubchannelToPointMapParameters>& params);

    //! Deconstructor
    virtual ~SubchannelToPointMap() { }

    /**
     * \brief Perform the map
     * \details  This performs the map of the output propertiy (Density or Temperature) to the given points.
     * \param f     Unused input
     * \param u     Vector containing the subchannel solution (may be a multivector containing other variables)
     * \param r     Vector to fill the densities.  The local size must match the number of points set
     *              by the call to setPoints.
     * \param a     Unused
     * \param b     Unused
     */
    void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
        AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0);

    virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
        return AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("Flow"));
    }


    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() { return d_outputVar; }


private:

    AMP_MPI d_comm;
    std::vector<double> d_point_x, d_point_y, d_point_z;
    AMP::LinearAlgebra::Variable::shared_ptr d_outputVar;

    // Create the subchannel grid for all processors
    size_t N_subchannels;
    std::vector<double> d_subchannel_x;     // x-coordinates of the center of the subchannel
    std::vector<double> d_subchannel_y;     // y-coordinates of the center of the subchannel
    std::vector<double> d_subchannel_z;     // z-coordinates of the z-faces of the subchannel
    void createGrid();

    AMP::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;
};


}
}

#endif

