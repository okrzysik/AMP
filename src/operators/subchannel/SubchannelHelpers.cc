#include "operators/subchannel/SubchannelHelpers.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "utils/Utilities.h"

#include <math.h>

namespace AMP {
namespace Operator {
namespace Subchannel {

// Compute the heat flux for the subchannel assuming a heat generation rate
std::vector<double> getHeatFluxGeneration( std::string heatShape, std::vector<double> z, double diameter, double Q_tot )
{
    for (size_t i=1; i<z.size(); i++)
        AMP_ASSERT(z[i]>z[i-1]);
    double height = z[z.size()-1];
    std::vector<double> dz(z.size()-1,0.0);
    for (size_t i=0; i<dz.size(); i++)
        dz[i] = z[i+1]-z[i];
    const double pi = 3.1415926535897932;
    std::vector<double> flux(dz.size(),0.0);
    if (heatShape == "Sinusoidal") {
        // sinusoidal
        for (size_t i=0; i<dz.size(); i++)
            flux[i] = Q_tot/(2.0*pi*diameter*dz[i]) * (cos(pi*z[i]/height) - cos(pi*z[i+1]/height));
    } else {
        AMP_ERROR("Heat shape '"+heatShape+" is invalid");
    }
    return flux;
}

// Compute the heat flux for the subchannel using the clad temperature
std::vector<double> getHeatFluxClad( std::vector<double> z, std::vector<AMP::Mesh::MeshElementID> face_ids, double diameter, 
    double channelDiam, double reynolds, double prandtl, double fraction, boost::shared_ptr<SubchannelPhysicsModel> subchannelPhysicsModel, 
    AMP::LinearAlgebra::Vector::const_shared_ptr flow, AMP::LinearAlgebra::Vector::const_shared_ptr clad_temp )
{
    for (size_t i=1; i<z.size(); i++)
        AMP_ASSERT(z[i]>z[i-1]);
    double height = z[z.size()-1];
    std::vector<double> dz(z.size()-1,0.0);
    for (size_t i=0; i<dz.size(); i++)
        dz[i] = z[i+1]-z[i];
    const double pi = 3.1415926535897932;
    AMP_ASSERT(face_ids.size()==z.size());
    AMP_ASSERT(flow!=NULL);
    AMP_ASSERT(clad_temp!=NULL);
    AMP::Discretization::DOFManager::shared_ptr flow_manager = flow->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr clad_manager = clad_temp->getDOFManager();
    const double h_scale = 1.0/Subchannel::scaleEnthalpy;   // Scale to change the input vector back to correct units
    const double P_scale = 1.0/Subchannel::scaleEnthalpy;   // Scale to change the input vector back to correct units 

    // Get the enthalapy, pressure, flow temperature, and clad temperature at the faces
    boost::shared_ptr<std::vector<double> >  h(new std::vector<double>(z.size(),0.0));
    boost::shared_ptr<std::vector<double> >  P(new std::vector<double>(z.size(),0.0));
    boost::shared_ptr<std::vector<double> >  Tf(new std::vector<double>(z.size(),0.0));
    boost::shared_ptr<std::vector<double> >  Tc(new std::vector<double>(z.size(),0.0));
    std::vector<size_t> flow_dofs(2), clad_dofs(1);
    for (size_t i=0; i<z.size(); i++) {
       flow_manager->getDOFs( face_ids[i], flow_dofs );
       clad_manager->getDOFs( face_ids[i], clad_dofs );
       AMP_ASSERT(flow_dofs.size()==2);
       AMP_ASSERT(clad_dofs.size()==1);
       (*h)[i] = h_scale*flow->getValueByGlobalID(flow_dofs[0]);
       (*P)[i] = P_scale*flow->getValueByGlobalID(flow_dofs[1]);
       (*Tc)[i] = clad_temp->getValueByGlobalID(clad_dofs[0]);
    }
    std::map<std::string, boost::shared_ptr<std::vector<double> > > temperatureArgMap;
    temperatureArgMap.insert(std::make_pair("enthalpy",h));
    temperatureArgMap.insert(std::make_pair("pressure",P));
    subchannelPhysicsModel->getProperty("Temperature", *Tf, temperatureArgMap);
    // Get the properties at cell centers
    size_t N = dz.size();
    boost::shared_ptr<std::vector<double> >  flowTemp(new std::vector<double>(N));
    boost::shared_ptr<std::vector<double> >  cladTemp(new std::vector<double>(N));
    boost::shared_ptr<std::vector<double> >  flowDens(new std::vector<double>(N));
    std::vector<double> specificVolume(z.size(),0.0);
    subchannelPhysicsModel->getProperty("SpecificVolume", specificVolume, temperatureArgMap);
    for (size_t i=0; i<N; i++) {
        (*flowTemp)[i] = 0.5*((*Tf)[i]+(*Tf)[i+1]);
        (*cladTemp)[i] = 0.5*((*Tc)[i]+(*Tc)[i+1]);
        (*flowDens)[i] = 0.5*(1./specificVolume[i]+1./specificVolume[+1]);
    }
    std::map<std::string, boost::shared_ptr<std::vector<double> > > convectiveHeatArgMap;
    convectiveHeatArgMap.insert(std::make_pair("temperature",cladTemp));
    convectiveHeatArgMap.insert(std::make_pair("density",flowDens));
    convectiveHeatArgMap.insert(std::make_pair("diameter",new std::vector<double>(N,channelDiam)));
    convectiveHeatArgMap.insert(std::make_pair("reynolds",new std::vector<double>(N,reynolds)));
    convectiveHeatArgMap.insert(std::make_pair("prandtl",new std::vector<double>(N,prandtl)));
    std::vector<double> heff(N); 
    subchannelPhysicsModel->getProperty("ConvectiveHeat", heff, convectiveHeatArgMap); 
    std::vector<double> flux(dz.size(),0.0);
    for (size_t i=0; i<N; i++) {
        double dT = (*cladTemp)[i] - (*flowTemp)[i];
        //flux[i] = heff[i]*dT*pi*diameter*fraction;
        flux[i] = heff[i]*dT*fraction;
    }
    return flux;
}



}
}
}


