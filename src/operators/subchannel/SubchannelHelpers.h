// This file containts some helper functions for subchannel operators
#ifndef included_AMP_SubchannelHelpers
#define included_AMP_SubchannelHelpers

#include "ampmesh/MeshElement.h"
#include "vectors/Vector.h"
#include "operators/subchannel/SubchannelPhysicsModel.h"

#include <string>
#include <vector>

namespace AMP {
namespace Operator {
namespace Subchannel {


/**
  * \brief Function to get the heat flux of the rod
  * \details  This function returns the heat flux of the rod (W/m^2) assuming a given generation rate
  * \param shape    The heat shape 
  * \param z        The axial positions of the faces (m)
  * \param diam     The diameter of the fuel rods (m)
  * \param Q_tot    The total heat generation rate (W)
  */
std::vector<double> getHeatFluxGeneration( std::string shape, std::vector<double> z, double diam, double Q_tot );


/**
  * \brief Function to get the heat flux of the rod
  * \details  This function returns the heat flux of the rod (W/m^2)
  * \param z            The axial positions of the faces
  * \param ids          Element ids of the faces of interest (only used for source="averageCladdingTemperature")
  * \param diam         The diameter of the fuel rods
  * \param channelDiam  The effective channel diameter
  * \param reynolds     The reynolds number
  * \param prandtl      The prandtl number
  * \param fraction     The fraction of a rod in the channel
  * \param subchannelPhysicsModel  The subchannel physics model
  * \param flow         The flow vector (h and P) (only used for source="averageCladdingTemperature")
  * \param clad_temp    The clad temperature mapped onto the faces (only used for source="averageCladdingTemperature")
  */
std::vector<double> getHeatFluxClad( std::vector<double> z, std::vector<AMP::Mesh::MeshElementID> face_ids, double diam, 
    double channelDiam, double reynolds, double prandtl, double fraction, boost::shared_ptr<SubchannelPhysicsModel> subchannelPhysicsModel, 
    AMP::LinearAlgebra::Vector::const_shared_ptr flow, AMP::LinearAlgebra::Vector::const_shared_ptr clad_temp );


}
}
}

#endif

