#ifndef included_AMP_SubchannelDensityToPointMap
#define included_AMP_SubchannelDensityToPointMap

#include "utils/AMP_MPI.h"
#include "ampmesh/Mesh.h"
#include "operators/Operator.h"
#include "operators/subchannel/SubchannelDensityToPointMapParameters.h"
#include "operators/subchannel/SubchannelPhysicsModel.h"

namespace AMP {
namespace Operator {


/**
 * \class SubchannelDensityToPointMap
 * \brief A class used to map subchannel density to points
 *
 * \details  This class provides routines for mapping the subchannel flow density 
 *     to a set of points provided.
 */
class SubchannelDensityToPointMap : public AMP::Operator::Operator
{
public :

    //! Default constructor
    SubchannelDensityToPointMap(const boost::shared_ptr<SubchannelDensityToPointMapParameters>& params);

    //! Deconstructor
    virtual ~SubchannelDensityToPointMap() { }

    /**
     * \brief Perform the map
     * \details  This performs the map of the density to the given points.
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


    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        return AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("Density"));
    }


private:

    AMP_MPI d_comm;
    std::vector<double> d_point_x, d_point_y, d_point_z;

    // Create the subchannel grid for all processors
    size_t N_subchannels;
    std::vector<double> d_subchannel_x;     // x-coordinates of the center of the subchannel
    std::vector<double> d_subchannel_y;     // y-coordinates of the center of the subchannel
    std::vector<double> d_subchannel_z;     // z-coordinates of the z-faces of the subchannel
    void createGrid();

    boost::shared_ptr<SubchannelPhysicsModel> d_subchannelPhysicsModel;
};


}
}

#endif

