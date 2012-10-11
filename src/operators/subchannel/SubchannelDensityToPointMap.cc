#include "operators/subchannel/SubchannelDensityToPointMap.h"
#include "utils/Utilities.h"
#include "ampmesh/StructuredMeshHelper.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "utils/ProfilerApp.h"


namespace AMP {
namespace Operator {



// Constructor
SubchannelDensityToPointMap::SubchannelDensityToPointMap(const boost::shared_ptr<SubchannelDensityToPointMapParameters>& params)
{
    // Copy the inputs
    AMP_ASSERT(params!=NULL);
    d_Mesh = params->d_Mesh;
    d_comm = params->d_comm;
    d_point_x = params->x;
    d_point_y = params->y;
    d_point_z = params->z;
    d_subchannelPhysicsModel = params->d_subchannelPhysicsModel;
    // Check the inputs
    AMP_INSIST(!d_comm.isNull(),"d_comm must not by NULL");
    AMP_INSIST(d_comm.anyReduce(d_Mesh!=NULL),"d_Mesh must be set on at least one processor");
    AMP_MPI mesh_comm(AMP_COMM_SELF);
    if ( d_Mesh!=NULL )
        mesh_comm = d_Mesh->getComm();
    AMP_INSIST(d_comm>=mesh_comm,"d_comm must be >= comm of the subchannel mesh");
    AMP_ASSERT(d_point_x.size()==d_point_y.size()&&d_point_x.size()==d_point_z.size());
    AMP_ASSERT(d_subchannelPhysicsModel!=NULL);
    // Get the coordinates of the subchannel mesh
    this->createGrid();
}


// Perform the map
void SubchannelDensityToPointMap::apply(AMP::LinearAlgebra::Vector::const_shared_ptr, AMP::LinearAlgebra::Vector::const_shared_ptr u,
    AMP::LinearAlgebra::Vector::shared_ptr r, const double a, const double b)
{
    PROFILE_START("apply");
    // Get the list of all subchannel face densities
    std::vector<double> density(N_subchannels*d_subchannel_z.size(),0.0);
    if ( d_Mesh!=NULL ) {
        AMP::LinearAlgebra::Vector::const_shared_ptr uInternal = subsetInputVector( u );
        AMP_ASSERT(uInternal!=NULL);
        AMP::Discretization::DOFManager::shared_ptr faceDOFManager = uInternal->getDOFManager(); 
        AMP::Mesh::MeshIterator it = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(d_Mesh,0);
        std::vector<size_t> dofs;
        const double h_scale = 1.0/Subchannel::scaleEnthalpy;                 // Scale to change the input vector back to correct units
        const double P_scale = 1.0/Subchannel::scalePressure;                 // Scale to change the input vector back to correct units
        for (size_t i=0; i<it.size(); i++) {
            // Get the density for the current face
            faceDOFManager->getDOFs( it->globalID(), dofs );
            AMP_ASSERT(dofs.size()==2);
            std::map<std::string, boost::shared_ptr<std::vector<double> > > subchannelArgMap;
            subchannelArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_scale*uInternal->getValueByGlobalID(dofs[0]))));
            subchannelArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,P_scale*uInternal->getValueByGlobalID(dofs[1]))));
            std::vector<double> specificVolume(1);
            d_subchannelPhysicsModel->getProperty("SpecificVolume", specificVolume, subchannelArgMap);
            double density_face = 1.0/specificVolume[0];
            // Add it to the subchannel density vector
            std::vector<double> center = it->centroid();
            size_t ix = AMP::Utilities::findfirst( d_subchannel_x, center[0]-1e-9 );
            size_t iy = AMP::Utilities::findfirst( d_subchannel_y, center[1]-1e-9 );
            size_t iz = AMP::Utilities::findfirst( d_subchannel_z, center[2]-1e-9 );
            AMP_ASSERT(AMP::Utilities::approx_equal(d_subchannel_x[ix],center[0]));
            AMP_ASSERT(AMP::Utilities::approx_equal(d_subchannel_y[iy],center[1]));
            AMP_ASSERT(AMP::Utilities::approx_equal(d_subchannel_z[iz],center[2]));
            size_t index = ix + iy*d_subchannel_x.size() + iz*N_subchannels;
            density[index] += density_face;
            ++it;
        }
    }
    d_comm.sumReduce<double>(&density[0],density.size());

    // Perform tri-linear interpolation to fill the density
    AMP::LinearAlgebra::VS_Comm commSelector( d_comm );
    AMP::LinearAlgebra::Vector::shared_ptr  densityVec = r->select(commSelector, u->getVariable()->getName());
    if ( densityVec!=NULL )
        densityVec = densityVec->subsetVectorForVariable(getOutputVariable());
    std::vector<double> localDensity(d_point_x.size(),0.0);
    if ( d_subchannel_x.size()>1 && d_subchannel_y.size()>1 ) {
        for (size_t i=0; i<d_point_x.size(); i++)
            localDensity[i] = AMP::Utilities::trilinear( d_subchannel_x, d_subchannel_y, 
                d_subchannel_z, density, d_point_x[i], d_point_y[i], d_point_z[i] );
    } else if ( d_subchannel_x.size()>1 ) {
        for (size_t i=0; i<d_point_x.size(); i++)
            localDensity[i] = AMP::Utilities::bilinear( d_subchannel_x, d_subchannel_z, 
                density, d_point_x[i], d_point_z[i] );
    } else if ( d_subchannel_y.size()>1 ) {
        for (size_t i=0; i<d_point_x.size(); i++)
            localDensity[i] = AMP::Utilities::bilinear( d_subchannel_y, d_subchannel_z, 
                density, d_point_y[i], d_point_z[i] );
    } else  {
        for (size_t i=0; i<d_point_x.size(); i++)
            localDensity[i] = AMP::Utilities::linear( d_subchannel_z, density, d_point_z[i] );
    }
    if ( d_point_x.size()>0 ) {
        AMP_ASSERT(densityVec!=NULL);
        AMP_ASSERT(densityVec->getLocalSize()==d_point_x.size());
        std::vector<size_t> dofs(d_point_x.size());
        for (size_t i=0; i<d_point_x.size(); i++)
            dofs[i] = i;
        densityVec->setValuesByLocalID( dofs.size(), &dofs[0], &localDensity[0] );
    }
    if ( densityVec )
        densityVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    PROFILE_STOP("apply");
}


// Create the subchannel grid for all processors
void SubchannelDensityToPointMap::createGrid()
{
    PROFILE_START("createGrid");
    std::set<double> x, y, z;
    if ( d_Mesh.get() != NULL ) {
        AMP::Mesh::MeshIterator it = d_Mesh->getIterator( AMP::Mesh::Vertex, 0 );
        for (size_t i=0; i<it.size(); i++) {
            std::vector<double> coord = it->coord();
            AMP_ASSERT(coord.size()==3);
            x.insert( coord[0] );
            y.insert( coord[1] );
            z.insert( coord[2] );
            ++it;
        }
    }
    d_comm.setGather(x);
    d_comm.setGather(y);
    d_comm.setGather(z);
    double last = 1e300;
    for (std::set<double>::iterator it=x.begin(); it!=x.end(); ++it) {
        if ( Utilities::approx_equal(last,*it,1e-12) )
            x.erase(it);
        else
            last = *it;
    }
    for (std::set<double>::iterator it=y.begin(); it!=y.end(); ++it) {
        if ( Utilities::approx_equal(last,*it,1e-12) )
            y.erase(it);
        else
            last = *it;
    }
    for (std::set<double>::iterator it=z.begin(); it!=z.end(); ++it) {
        if ( Utilities::approx_equal(last,*it,1e-12) )
            z.erase(it);
        else
            last = *it;
    }
    std::vector<double> x_grid(x.begin(),x.end());
    d_subchannel_x = std::vector<double>(x_grid.size()-1,0.0);
    for (size_t i=0; i<d_subchannel_x.size(); i++)
        d_subchannel_x[i] = 0.5*(x_grid[i]+x_grid[i+1]);
    std::vector<double> y_grid(y.begin(),y.end());
    d_subchannel_y = std::vector<double>(x_grid.size()-1,0.0);
    for (size_t i=0; i<d_subchannel_y.size(); i++)
        d_subchannel_y[i] = 0.5*(y_grid[i]+y_grid[i+1]);
    d_subchannel_z = std::vector<double>(z.begin(),z.end());

    size_t Nx = d_subchannel_x.size();
    size_t Ny = d_subchannel_y.size();
    size_t Nz = d_subchannel_z.size()-1;
    if ( d_Mesh.get() != NULL ) 
        AMP_ASSERT(Nx*Ny*Nz==d_Mesh->numGlobalElements(AMP::Mesh::Volume));
    N_subchannels = Nx*Ny;
    PROFILE_STOP("createGrid");
}


}
}

