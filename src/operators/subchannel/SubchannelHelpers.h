// This file containts some helper functions for subchannel operators
#ifndef included_AMP_SubchannelHelpers
#define included_AMP_SubchannelHelpers

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "vectors/Vector.h"
#include "operators/subchannel/SubchannelPhysicsModel.h"

#include <string>
#include <vector>

namespace AMP {
namespace Operator {
namespace Subchannel {


/**
  * \brief Function to get some basic properties for each subchannel based on the mesh
  * \details  This function returns some basic properties including the hydraulic diameter, subchannel area, etc
  *   from the mesh
  * \param[in] subchannel   The subchannel mesh
  * \param[in] clad_x       The x-coordinates of the clad
  * \param[in] clad_y       The y-coordinates of the clad
  * \param[in] clad_d       The diameters of the clad
  * \param[out] x           The x-coordinates of the subchannel boundaries (Nx+1)
  * \param[out] y           The y-coordinates of the subchannel boundaries (Ny+1)
  * \param[out] area        The flow area of the subchannels (Nx x Ny)
  * \param[out] diam        The hydraulic diameter of the subchannels (Nx x Ny)
  * \param[out] perimeter   The wetted perimeter of the subchannels (Nx x Ny)
  * \param[out] rod_diameter  The average rod diameter for each subchannel (Nx x Ny)
  * \param[out] channel_fraction  The fraction of the rod in each subchannel (Nx x Ny)
  */
void getSubchannelProperties( AMP::Mesh::Mesh::shared_ptr subchannel, const std::vector<double>& clad_x,
     const std::vector<double>& clad_y, const std::vector<double>& clad_d, std::vector<double>& x,
     std::vector<double>& y, std::vector<double>& area, std::vector<double>& diam, std::vector<double>& perimeter, 
     std::vector<double>& rod_diameter, std::vector<double>& channel_fraction );


/**
  * \brief Function to get the coordinates of the clad
  * \details  This function returns the coordinates and diameter of the clad meshes across all processors
  * \param[in] comm         Communicator over which we want to gather the results (must be >= clad mesh comm)
  * \param[in] clad         Multimesh containing the clad meshes
  * \param[out] x           The x-coordinates of the clad
  * \param[out] y           The y-coordinates of the clad
  * \param[out] diam        The diameters of the clad
  */
void getCladProperties( AMP::AMP_MPI comm, AMP::Mesh::Mesh::shared_ptr clad, std::vector<double>& x,
     std::vector<double>& y, std::vector<double>& diam );


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
  * \param face_ids     Element ids of the faces of interest (only used for source="averageCladdingTemperature")
  * \param channelDiam  The effective channel diameter
  * \param reynolds     The reynolds number
  * \param prandtl      The prandtl number
  * \param fraction     The fraction of a rod in the channel
  * \param subchannelPhysicsModel  The subchannel physics model
  * \param flow         The flow vector (h and P) (only used for source="averageCladdingTemperature")
  * \param clad_temp    The clad temperature mapped onto the faces (only used for source="averageCladdingTemperature")
  */
std::vector<double> getHeatFluxClad( std::vector<double> z, std::vector<AMP::Mesh::MeshElementID> face_ids,
    double channelDiam, double reynolds, double prandtl, double fraction, boost::shared_ptr<SubchannelPhysicsModel> subchannelPhysicsModel, 
    AMP::LinearAlgebra::Vector::const_shared_ptr flow, AMP::LinearAlgebra::Vector::const_shared_ptr clad_temp );


}
}
}

#endif

