#include "operators/subchannel/SubchannelHelpers.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "ampmesh/StructuredMeshHelper.h"
#include "utils/Utilities.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include <math.h>

namespace AMP {
namespace Operator {
namespace Subchannel {



// Compute basic properties from the subchannel mesh
void getSubchannelProperties( AMP::Mesh::Mesh::shared_ptr subchannel, const std::vector<double>& clad_x,
     const std::vector<double>& clad_y, const std::vector<double>& clad_d, std::vector<double>& x,
     std::vector<double>& y, std::vector<double>& area, std::vector<double>& fric_diam, std::vector<double>& heat_diam,
     std::vector<double>& rod_diameter, std::vector<double>& channel_fraction )
{
    const double pi = 3.1415926535897932;
    AMP_MPI comm = subchannel->getComm();
    // First get the x-y-z coordinates of the subchannel mesh
    std::vector<double> z;
    AMP::Mesh::StructuredMeshHelper::getXYZCoordinates( subchannel, x, y, z );
    AMP_ASSERT(x.size()>=2&&y.size()>=2);
    size_t Nx = x.size()-1;
    size_t Ny = y.size()-1;
    // Check the clad properties
    AMP_ASSERT(clad_x.size()==clad_y.size()&&clad_x.size()==clad_d.size());
    AMP_ASSERT(clad_x.size()==comm.maxReduce(clad_x.size()));
    AMP_ASSERT(clad_x.size()>0);
    // For each subchannel, get the area, wetted perimeter, and the hydraulic diameter
    size_t N_subchannels = (x.size()-1)*(y.size()-1);
    area = std::vector<double>(N_subchannels,0.0);
    fric_diam = std::vector<double>(N_subchannels,0.0);
    heat_diam = std::vector<double>(N_subchannels,0.0);
    rod_diameter = std::vector<double>(N_subchannels,0.0);
    channel_fraction = std::vector<double>(N_subchannels,0.0);
    std::vector<double> perimeter1(N_subchannels,0.0);    // Clad only perimeter
    std::vector<double> perimeter2(N_subchannels,0.0);    // Clad + water perimeter
    // Get the area of the subchannel without the clad
    for (size_t i=0; i<N_subchannels; i++) {
        size_t ix = i%(x.size()-1);
        size_t iy = i/(x.size()-1);
        area[i] = (x[ix+1]-x[ix])*(y[iy+1]-y[iy]);
        perimeter2[i] = 2*(x[ix+1]-x[ix]) + 2*(y[iy+1]-y[iy]);
    }
    // Add the area and perimeter corrections of the clad (assuming no clads overlap)
    const double TOL=1e-12;
    for (size_t k=0; k<clad_x.size(); k++) {
        double xc = clad_x[k];
        double yc = clad_y[k];
        double dc = clad_d[k];
        size_t index_x = AMP::Utilities::findfirst(x,xc-TOL);
        size_t index_y = AMP::Utilities::findfirst(y,yc-TOL);
        if ( index_x==x.size() ) { index_x--; }
        if ( index_y==y.size() ) { index_y--; }
        if ( fabs(x[index_x]-xc)<=TOL && fabs(y[index_y]-yc)<=TOL ) {
            // The clad is located at the subchannel boundaries
            double dA = 0.25*pi*0.25*dc*dc;
            double dP1 = 0.25*pi*dc;
            double dP2 = (1.0-0.25*pi)*dc;
            size_t i[4];
            for (int j=0; j<4; j++)
                i[j] = static_cast<size_t>(-1);
            if ( index_x>0 && index_y>0 )
                i[0] = index_x-1 + (index_y-1)*Nx;
            if ( index_x<Nx && index_y>0 )
                i[1] = index_x   + (index_y-1)*Nx;
            if ( index_x>0 && index_y<Ny )
                i[2] = index_x-1 + index_y*Nx;
            if ( index_x<Nx && index_y<Ny )
                i[3] = index_x   + index_y*Nx;
            for (int j=0; j<4; j++) {
                if ( i[j]==static_cast<size_t>(-1) )
                    continue;
                area[i[j]] -= dA;
                perimeter1[i[j]] += dP1;
                perimeter2[i[j]] -= dP2;
                double ratio = 1.0/(channel_fraction[i[j]]+1.0);
                channel_fraction[i[j]] += 0.25;
                rod_diameter[i[j]] = (1.0-ratio)*rod_diameter[i[j]] + ratio*dc;
            }
        } else {
            if ( index_x==Nx ) { index_x--; }
            if ( index_y==Nx ) { index_y--; }
            if ( (xc-0.5*dc)>x[index_x] && (xc+0.5*dc)<x[index_x+1] && 
                 (yc-0.5*dc)>y[index_y] && (yc+0.5*dc)<y[index_y+1] ) {
                // The clad inside the subchannel
                size_t i = index_x+index_y*Nx;
                area[i] -= 0.25*pi*dc*dc;
                perimeter1[i] += pi*dc;
                perimeter2[i] += pi*dc;
                double R = 1.0/(channel_fraction[i]+1.0);
                channel_fraction[i] += 1.0;
                rod_diameter[i] = (1.0-R)*rod_diameter[i] + R*dc;
            } else {
                AMP_ERROR("General case not handled yet\n");
            }
        }
    }
    // Compute the hydraulic diameter
    for (size_t i=0; i<N_subchannels; i++) {
        fric_diam[i] = 4.0*area[i]/perimeter1[i];
        heat_diam[i] = 4.0*area[i]/perimeter2[i];
    }
}


// Function to get the clad properties
void getCladProperties( AMP::AMP_MPI comm, AMP::Mesh::Mesh::shared_ptr clad, std::vector<double>& x,
     std::vector<double>& y, std::vector<double>& diam )
{
    // Get the center of each local clad
    std::set< AMP::Utilities::triplet<double,double,double> >  center;
    if ( clad!=NULL ) {
        AMP_ASSERT(clad->getComm()<=comm);
        std::vector<AMP::Mesh::MeshID> ids = clad->getLocalBaseMeshIDs();
        for (size_t i=0; i<ids.size(); i++) {
            AMP::Mesh::Mesh::shared_ptr mesh = clad->Subset(ids[i]);
            std::vector<double> box = mesh->getBoundingBox();
            AMP::Utilities::triplet<double,double,double> tmp;
            tmp.first = 0.5*(box[0]+box[1]);
            tmp.second = 0.5*(box[2]+box[3]);
            tmp.third = std::max( std::max(tmp.first-box[0],box[1]-tmp.first),
                                  std::max(tmp.second-box[2],box[3]-tmp.second) );
            center.insert( tmp );
        }
    }
    // Get the global set and check that there are no duplicates
    comm.setGather(center);
    std::vector< AMP::Utilities::triplet<double,double,double> >  center2(center.begin(),center.end());
    for (size_t i=0; i<center2.size(); i++) {
        for (size_t j=i+1; j<center2.size(); j++) {
            if ( AMP::Utilities::approx_equal(center2[i].first,center2[j].first) &&
                 AMP::Utilities::approx_equal(center2[i].second,center2[j].second) ) {
                AMP_ERROR("Duplicate clads detected");
            }
        }
    }
    x.resize(center2.size());
    y.resize(center2.size());
    diam.resize(center2.size());
    for (size_t i=0; i<center2.size(); i++) {
        x[i] = center2[i].first;
        y[i] = center2[i].second;
        diam[i] = 2.0*center2[i].third;
    }
}


// Compute the heat flux for the subchannel assuming a heat generation rate
std::vector<double> getHeatFluxGeneration( std::string heatShape, std::vector<double> z, double diameter, double Q_tot )
{
    for (size_t i=1; i<z.size(); i++)
        AMP_ASSERT(z[i]>z[i-1]);
    double height = z.back()-z.front();
    std::vector<double> dz(z.size()-1,0.0);
    for (size_t i=0; i<dz.size(); i++)
        dz[i] = z[i+1]-z[i];
    const double pi = 3.1415926535897932;
    std::vector<double> flux(dz.size(),0.0);
    if (heatShape == "Flat") {
        // sinusoidal
        for (size_t i=0; i<dz.size(); i++)
            flux[i] = Q_tot/(pi*diameter*height);
    } else if (heatShape == "Sinusoidal") {
        // sinusoidal
        for (size_t i=0; i<dz.size(); i++)
            flux[i] = Q_tot/(2.0*pi*diameter*dz[i]) * (cos(pi*(z[i]-z[0])/height) - cos(pi*(z[i+1]-z[0])/height));
    } else {
        AMP_ERROR("Heat shape '"+heatShape+" is invalid");
    }
    return flux;
}

// Compute the heat flux for the subchannel using the clad temperature
std::vector<double> getHeatFluxClad( std::vector<double> z, std::vector<AMP::Mesh::MeshElementID> face_ids,
    double channelDiam, double reynolds, double prandtl, double fraction, boost::shared_ptr<SubchannelPhysicsModel> subchannelPhysicsModel, 
    AMP::LinearAlgebra::Vector::const_shared_ptr flow, AMP::LinearAlgebra::Vector::const_shared_ptr clad_temp )
{
    for (size_t i=1; i<z.size(); i++)
        AMP_ASSERT(z[i]>z[i-1]);
    std::vector<double> dz(z.size()-1,0.0);
    for (size_t i=0; i<dz.size(); i++)
        dz[i] = z[i+1]-z[i];
    //const double pi = 3.1415926535897932;
    AMP_ASSERT(face_ids.size()==z.size());
    AMP_ASSERT(flow!=NULL);
    AMP_ASSERT(clad_temp!=NULL);
    AMP::Discretization::DOFManager::shared_ptr flow_manager = flow->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr clad_manager = clad_temp->getDOFManager();
    const double h_scale = 1.0/Subchannel::scaleEnthalpy;   // Scale to change the input vector back to correct units
    const double P_scale = 1.0/Subchannel::scalePressure;   // Scale to change the input vector back to correct units 

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
        flux[i] = heff[i]*dT*fraction;
    }
    return flux;
}


// Function to compute the hydraulic diameter of the subchannels on the clad surface
AMP::LinearAlgebra::Vector::shared_ptr  getCladHydraulicDiameter( AMP::Mesh::Mesh::shared_ptr clad, 
    AMP::Mesh::Mesh::shared_ptr subchannel, AMP::AMP_MPI comm )
{
    if ( clad.get() != NULL )
        AMP_ASSERT(clad->getComm()<=comm);
    if ( subchannel.get() != NULL )
        AMP_ASSERT(subchannel->getComm()<=comm);
    // Get the clad properties
    std::vector<double> clad_x, clad_y, clad_d;
    getCladProperties( comm, clad, clad_x, clad_y, clad_d );
    AMP::Mesh::Mesh::shared_ptr clad_surface;
    if ( clad.get()!=NULL )
        clad_surface = clad->Subset( clad->getBoundaryIDIterator(AMP::Mesh::Face,4,1) );
    // Get the subchannel properties
    size_t N_subchannels;
    std::vector<double> x, y, heat_diam;
    int root = -1;
    if ( subchannel.get() != NULL ) {
        std::vector<double> area, fric_diam, rod_diameter, channel_fraction;
        getSubchannelProperties( subchannel, clad_x, clad_y, clad_d, 
            x, y, area, fric_diam, heat_diam, rod_diameter, channel_fraction );
        N_subchannels = x.size();
        root = comm.getRank();
    }
    root = comm.maxReduce(root);
    N_subchannels = comm.bcast(N_subchannels,root);
    if ( subchannel.get() == NULL ) {
        x.resize(N_subchannels);
        y.resize(N_subchannels);
        heat_diam.resize(N_subchannels);
    }
    comm.bcast(&x[0],N_subchannels,root);
    comm.bcast(&y[0],N_subchannels,root);
    comm.bcast(&heat_diam[0],N_subchannels,root);
    // Return if we are not on the clad surface
    if ( clad_surface.get()==NULL )
        return AMP::LinearAlgebra::Vector::shared_ptr();
    // Create and initialize the vector
    AMP::Discretization::DOFManager::shared_ptr DOF = 
        AMP::Discretization::simpleDOFManager::create(clad_surface,AMP::Mesh::Vertex,1,1,true);
    AMP::LinearAlgebra::Variable::shared_ptr variable( new AMP::LinearAlgebra::Variable("ChannelDiameter") );
    AMP::LinearAlgebra::Vector::shared_ptr diameter = AMP::LinearAlgebra::createVector( DOF, variable );
    diameter->zero();
    AMP::Mesh::MeshIterator it = clad_surface->getIterator(AMP::Mesh::Vertex);
    std::vector<size_t> dofs(1);
    size_t Nx = x.size()-1;
    size_t Ny = y.size()-1;
    for (size_t i=0; i<it.size(); i++) {
        std::vector<double> pos = it->coord();
        DOF->getDOFs(it->globalID(),dofs);
        AMP_ASSERT(dofs.size()==1);
        size_t ix = AMP::Utilities::findfirst(x,pos[0]);
        size_t iy = AMP::Utilities::findfirst(x,pos[0]);
        if ( ix==0 ) { ix = 1; }
        if ( iy==0 ) { iy = 1; }
        ix--;
        iy--;
        if ( ix==Nx ) { ix = Nx-1; }
        if ( iy==Nx ) { iy = Ny-1; }
        diameter->setValueByGlobalID( dofs[0], heat_diam[ix+iy*Nx] );
        ++it;
    }
    diameter->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
    return diameter;
}


}
}
}


