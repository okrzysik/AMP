
#include "SubchannelFourEqNonlinearOperator.h"
#include "SubchannelOperatorParameters.h"
#include "ampmesh/StructuredMeshHelper.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

#include <string>


namespace AMP {
namespace Operator {

// reset
void SubchannelFourEqNonlinearOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
{
      boost::shared_ptr<SubchannelOperatorParameters> myparams = 
        boost::dynamic_pointer_cast<SubchannelOperatorParameters>(params);

      AMP_INSIST( ((myparams.get()) != NULL), "NULL parameters" );
      AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );

      d_Pout = getDoubleParameter(myparams,"Exit_Pressure",15.5132e6);
      d_Tin  = getDoubleParameter(myparams,"Inlet_Temperature",569.26);  
      d_min  = getDoubleParameter(myparams,"Inlet_Axial_Flow_Rate",0.3522);  
      d_win  = getDoubleParameter(myparams,"Inlet_Lateral_Flow_Rate",0.0);  
      d_area = (myparams->d_db)->getDoubleArray("SubchannelArea");
      d_gamma     = getDoubleParameter(myparams,"Fission_Heating_Coefficient",0.0);  
      d_theta     = getDoubleParameter(myparams,"Channel_Angle",0.0);  
      d_friction  = getDoubleParameter(myparams,"Friction_Factor",0.1);  
      d_pitch     = getDoubleParameter(myparams,"Lattice_Pitch",0.0128016);  
      d_diameter  = getDoubleParameter(myparams,"Rod_Diameter",0.0097028);  
      d_K    = getDoubleParameter(myparams,"Form_Loss_Coefficient",0.2);  
      d_source = getStringParameter(myparams,"Heat_Source_Type","totalHeatGeneration");

      // get additional parameters based on heat source type
      if (d_source == "totalHeatGeneration") {
          d_Q    = getDoubleParameter(myparams,"Rod_Power",66.81e3);  
          d_heatShape = getStringParameter(myparams,"Heat_Shape","Sinusoidal");
      }

      d_subchannelPhysicsModel = myparams->d_subchannelPhysicsModel;  
}

// function used in reset to get double parameter or set default if missing
double SubchannelFourEqNonlinearOperator::getDoubleParameter(	boost::shared_ptr<SubchannelOperatorParameters> myparams,
								std::string paramString,
                                                                double defaultValue)
{
    bool keyExists = (myparams->d_db)->keyExists(paramString);
    if (keyExists) {
       return (myparams->d_db)->getDouble(paramString);
    } else {
       AMP::pout << "Key '"+paramString+"' was not provided. Using default value: " << defaultValue << "\n";
       return defaultValue;
    }
}

// function used in reset to get string parameter or set default if missing
std::string SubchannelFourEqNonlinearOperator::getStringParameter(	boost::shared_ptr<SubchannelOperatorParameters> myparams,
									std::string paramString,
                                                                	std::string defaultValue)
{
    bool keyExists = (myparams->d_db)->keyExists(paramString);
    if (keyExists) {
       return (myparams->d_db)->getString(paramString);
    } else {
       AMP::pout << "Key '"+paramString+"' was not provided. Using default value: " << defaultValue << "\n";
       return defaultValue;
    }
}

// function used to get all lateral gaps
std::map<std::vector<double>,AMP::Mesh::MeshElement> SubchannelFourEqNonlinearOperator::getLateralFaces(AMP::Mesh::Mesh::shared_ptr mesh)
{
   // map of lateral gaps to their centroids
   std::map<std::vector<double>,AMP::Mesh::MeshElement> lateralFaceMap;
   // get iterator over all faces of mesh
   AMP::Mesh::MeshIterator face = mesh->getIterator(AMP::Mesh::Face,0);
   // loop over faces
   for (; face != face.end(); face++) {
      // check that face is vertical
      // ---------------------------
      // get centroid of current face
      std::vector<double> faceCentroid = face->centroid();
      // get vertices of current face
      std::vector<AMP::Mesh::MeshElement> vertices = face->getElements(AMP::Mesh::Vertex);

      bool perpindicular_to_x = true; // is the current face perpindicular to x-axis?
      bool perpindicular_to_y = true; // is the current face perpindicular to y-axis?
      // loop over vertices of current face
      for (size_t j=0; j<vertices.size(); ++j) {
        // get coordinates of current vertex
        std::vector<double> vertexCoord = vertices[j].coord();
        // if any vertex does not have the same x-coordinate as the face centroid,
        if ( !AMP::Utilities::approx_equal(vertexCoord[0],faceCentroid[0], 1.0e-6) )
          // then the face is not perpindicular to x-axis
          perpindicular_to_x = false;
        // if any vertex does not have the same y-coordinate as the face centroid,
        if ( !AMP::Utilities::approx_equal(vertexCoord[1],faceCentroid[1], 1.0e-6) )
          // then the face is not perpindicular to y-axis
          perpindicular_to_y = false;
      }
      // check that face is in the interior of the mesh; it must have two adjacent cells
      // -------------------------------------------------------------------------------
      // if the face is vertical
      if (perpindicular_to_x || perpindicular_to_y) {
         // if the face has more than 1 adjacent cell
         if ((mesh->getElementParents(*face,AMP::Mesh::Volume)).size() > 1) {
            // insert face into map with centroid
            lateralFaceMap.insert(std::pair<std::vector<double>,AMP::Mesh::MeshElement>(faceCentroid,*face));
         }
      }
   }// end loop over faces
   return lateralFaceMap;
}

// function used to get all of the unique x,y,z points in subchannel mesh
void SubchannelFourEqNonlinearOperator::fillSubchannelGrid(AMP::Mesh::Mesh::shared_ptr mesh)
{
    // Create the grid for all processors
    std::set<double> x, y, z;
    if ( mesh.get() != NULL ) {
        AMP::Mesh::MeshIterator vertex = mesh->getIterator( AMP::Mesh::Vertex, 0 );
        // for all vertices in mesh
        for (size_t i=0; i<vertex.size(); i++) {
            std::vector<double> coord = vertex->coord();
            AMP_ASSERT(coord.size()==3);
            // insert x,y,z points into sets, even if duplicate
            x.insert( coord[0] );
            y.insert( coord[1] );
            z.insert( coord[2] );
            ++vertex;
        }
    }
    d_Mesh->getComm().setGather(x);
    d_Mesh->getComm().setGather(y);
    d_Mesh->getComm().setGather(z);
    double last = 1.0e300; // arbitary large number
    // erase duplicate x points
    for (std::set<double>::iterator it=x.begin(); it!=x.end(); ++it) {
        if ( Utilities::approx_equal(last,*it,1e-12) )
            x.erase(it);
        else
            last = *it;
    }
    // erase duplicate y points
    for (std::set<double>::iterator it=y.begin(); it!=y.end(); ++it) {
        if ( Utilities::approx_equal(last,*it,1e-12) )
            y.erase(it);
        else
            last = *it;
    }
    // erase duplicate z points
    for (std::set<double>::iterator it=z.begin(); it!=z.end(); ++it) {
        if ( Utilities::approx_equal(last,*it,1e-12) )
            z.erase(it);
        else
            last = *it;
    }
    d_x = std::vector<double>(x.begin(),x.end());
    d_y = std::vector<double>(y.begin(),y.end());
    d_z = std::vector<double>(z.begin(),z.end());
    size_t Nx = d_x.size()-1; // number of mesh divisions along x-axis
    size_t Ny = d_y.size()-1; // number of mesh divisions along y-axis
    size_t Nz = d_z.size()-1; // number of mesh divisions along z-axis
    if ( mesh.get() != NULL ) 
        // check that computed number of elements matches that found by numGlobalElements()
        AMP_ASSERT(Nx*Ny*Nz==mesh->numGlobalElements(AMP::Mesh::Volume));
    // compute number of subchannels
    d_numSubchannels = Nx*Ny;
}

// function to get a unique index for a subchannel based on its x,y coordinates
int SubchannelFourEqNonlinearOperator::getSubchannelIndex( double x, double y )
{
    // get index of first entry in subchannel x mesh >= x
    size_t ix = Utilities::findfirst(d_x,x);
    // get index of first entry in subchannel y mesh >= y
    size_t iy = Utilities::findfirst(d_y,y);
    // check that indices are valid
    if ( ix>0 && ix<d_x.size() && iy>0 && iy<d_y.size() ) {
        return (ix-1)+(iy-1)*(d_x.size()-1);
    } else {
        AMP_ERROR("Invalid indices found for getSubchannelIndex()");
    }
    return 0;
}

// apply
void SubchannelFourEqNonlinearOperator :: apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
    AMP::LinearAlgebra::Vector::shared_ptr r, const double a, const double b)
{

      // ensure that solution and residual vectors aren't NULL
      AMP_INSIST( ((r.get()) != NULL), "NULL Residual Vector" );
      AMP_INSIST( ((u.get()) != NULL), "NULL Solution Vector" );

      // calculate extra parameters
      const double pi = 4.0*atan(1.0); // pi
      const double g = 9.805;          // acceleration due to gravity [m/s2]
      // assuming square pitch
      const double perimeter = 4.0*(d_pitch-d_diameter) + pi*d_diameter;    // wetted perimeter
      const double area = std::pow(d_pitch,2) - pi*std::pow(d_diameter,2)/4.0; // flow area
      const double D = 4.0*area/perimeter;                                     // hydraulic diameter

      // Subset the vectors
      AMP::LinearAlgebra::Vector::const_shared_ptr inputVec  = subsetInputVector( u );
      AMP::LinearAlgebra::Vector::shared_ptr       outputVec = subsetOutputVector( r );

      AMP::Discretization::DOFManager::shared_ptr dof_manager = inputVec->getDOFManager();

      // get all of the unique x,y,z points in subchannel mesh
      fillSubchannelGrid(d_Mesh);

      // get map of all of the lateral faces to their centroids
      std::map<std::vector<double>,AMP::Mesh::MeshElement> lateralFaceMap = getLateralFaces(d_Mesh);

      // compute height of subchannels
      std::vector<double> box = d_Mesh->getBoundingBox();
      const double height = box[5] - box[4];

      AMP::Mesh::MeshIterator cell = d_Mesh->getIterator(AMP::Mesh::Volume, 0); // iterator for cells of mesh

      std::vector<std::vector<AMP::Mesh::MeshElement> > d_elem(d_numSubchannels); // array of array of elements for each subchannel
      std::vector<bool> d_ownSubChannel(d_numSubchannels);

      // for each cell,
      for( ; cell != cell.end(); ++cell) {
        std::vector<double> center = cell->centroid();
        // get the index of the subchannel
        int index = getSubchannelIndex( center[0], center[1] );
        if ( index>=0 ){
          d_ownSubChannel[index] = true;
          // put cell into array of cells for that subchannel
          d_elem[index].push_back( *cell );
        }
      }// end for cell

      // for each subchannel,
      for(int isub =0; isub<d_numSubchannels; ++isub){
        if(d_ownSubChannel[isub]){
          // extract subchannel cells from d_elem[isub]
          boost::shared_ptr<std::vector<AMP::Mesh::MeshElement> > subchannelElements( new std::vector<AMP::Mesh::MeshElement>() );
          subchannelElements->reserve(d_numSubchannels);
          for(size_t ielem=0; ielem<d_elem[isub].size(); ++ielem){
            subchannelElements->push_back(d_elem[isub][ielem]);
          }
          AMP::Mesh::MeshIterator     localSubchannelCell = AMP::Mesh::MultiVectorIterator( subchannelElements ); // iterator over elements of current subchannel
          //AMP::Mesh::Mesh::shared_ptr localSubchannelMesh = d_Mesh->Subset( localSubchannelCell  ); // subset the mesh for the current subchannel
          // get subchannel index
          std::vector<double> subchannelCentroid = localSubchannelCell->centroid();
          size_t currentSubchannelIndex = getSubchannelIndex(subchannelCentroid[0],subchannelCentroid[1]);

          // loop over cells of current subchannel
          for (; localSubchannelCell != localSubchannelCell.end(); ++localSubchannelCell) {
             // get upper and lower axial faces of current cell
             // -----------------------------------------------
             std::vector<double> cellCentroid = localSubchannelCell->centroid();
             // get axial faces of current cell
             std::vector<AMP::Mesh::MeshElement> cellFaces = localSubchannelCell->getElements(AMP::Mesh::Face);
             AMP::Mesh::MeshElement plusFace; // upper axial face for current cell
             AMP::Mesh::MeshElement minusFace; // lower axial face for current cell
             // loop over faces of current cell
             for (std::vector<AMP::Mesh::MeshElement>::iterator face = cellFaces.begin(); face != cellFaces.end(); ++face) {
                std::vector<double> faceCentroid = face->centroid();
                // if z-coordinates of centroids of the cell and face are not equal,
                if (!AMP::Utilities::approx_equal(faceCentroid[2],cellCentroid[2],1.0e-6)) {
                   // if face is above cell centroid,
                   if (faceCentroid[2] > cellCentroid[2])
                      // face is the upper axial face
                      plusFace = *face;
                   else
                      // face is the lower axial face
                      minusFace = *face;
                }
             }
             std::vector<double> plusFaceCentroid = plusFace.centroid();
             std::vector<double> minusFaceCentroid = minusFace.centroid();

             // extract unknowns from solution vector
             // -------------------------------------
             std::vector<size_t> plusDofs;
             std::vector<size_t> minusDofs;
             dof_manager->getDOFs(plusFace.globalID(),plusDofs);
             dof_manager->getDOFs(minusFace.globalID(),minusDofs);
             double m_plus = inputVec->getValueByGlobalID(plusDofs[0]);
             double h_plus = inputVec->getValueByGlobalID(plusDofs[1]);
             double p_plus = inputVec->getValueByGlobalID(plusDofs[2]);
             double m_minus = inputVec->getValueByGlobalID(minusDofs[0]);
             double h_minus = inputVec->getValueByGlobalID(minusDofs[1]);
             double p_minus = inputVec->getValueByGlobalID(minusDofs[2]);

             // compute additional quantities
             // -----------------------------
             double m_mid = 1.0/2.0*(m_plus + m_minus);
             double p_mid = 1.0/2.0*(p_plus + p_minus);

             // evaluate specific volume at upper face
             std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_plus;
             volumeArgMap_plus.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_plus)));
             volumeArgMap_plus.insert(std::make_pair("pressure",new std::vector<double>(1,p_plus)));
             std::vector<double> volumeResult_plus(1);
             d_subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_plus,volumeArgMap_plus); 
             double vol_plus = volumeResult_plus[0];

             // evaluate specific volume at lower face
             std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_minus;
             volumeArgMap_minus.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_minus)));
             volumeArgMap_minus.insert(std::make_pair("pressure",new std::vector<double>(1,p_minus)));
             std::vector<double> volumeResult_minus(1);
             d_subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_minus,volumeArgMap_minus); 
             double vol_minus = volumeResult_minus[0];

             // determine axial donor quantities
             double h_axialDonor;
             double vol_axialDonor;
             if (m_mid >= 0.0) {
                vol_axialDonor = vol_minus;
                h_axialDonor = h_minus;
             } else {
                vol_axialDonor = vol_plus;
                h_axialDonor = h_minus;
             }

             double rho_mid = 1.0/vol_axialDonor;

             double u_plus  = m_plus*vol_plus/area;
             double u_minus = m_minus*vol_minus/area;
             double u_mid = m_mid*vol_axialDonor/area;

             // evaluate temperature for cell
             std::map<std::string, boost::shared_ptr<std::vector<double> > > temperatureArgMap_mid;
             temperatureArgMap_mid.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_axialDonor)));
             temperatureArgMap_mid.insert(std::make_pair("pressure",new std::vector<double>(1,p_mid)));
             std::vector<double> temperatureResult_mid(1);
             d_subchannelPhysicsModel->getProperty("Temperature",temperatureResult_mid,temperatureArgMap_mid); 
             double T_mid = temperatureResult_mid[0];

             // evaluate conductivity for cell
             std::map<std::string, boost::shared_ptr<std::vector<double> > > conductivityArgMap_mid;
             conductivityArgMap_mid.insert(std::make_pair("conductivity",new std::vector<double>(1,T_mid)));
             conductivityArgMap_mid.insert(std::make_pair("density",new std::vector<double>(1,rho_mid)));
             std::vector<double> conductivityResult_mid(1);
             d_subchannelPhysicsModel->getProperty("ThermalConductivity",conductivityResult_mid,conductivityArgMap_mid); 
             double k_mid = conductivityResult_mid[0];

             // compute element height
             double dz = plusFaceCentroid[2] - minusFaceCentroid[2];

             // initialize sum terms
             double mass_crossflow_sum = 0.0;
             double energy_crossflow_sum = 0.0;
             double energy_heatflux_sum = 0.0;
             double energy_turbulence_sum = 0.0;
             double energy_conduction_sum = 0.0;
             double energy_direct_heating_sum = 0.0;
             double axial_crossflow_sum = 0.0;
             double axial_turbulence_sum = 0.0;

             // loop over gap faces
             for (std::vector<AMP::Mesh::MeshElement>::iterator face = cellFaces.begin(); face != cellFaces.end(); ++face) {
                std::vector<double> faceCentroid = face->centroid();
                std::map<std::vector<double>,AMP::Mesh::MeshElement>::iterator lateralFaceIterator = lateralFaceMap.find(faceCentroid);
                if (lateralFaceIterator != lateralFaceMap.end()) {
                   // get face
                   AMP::Mesh::MeshElement lateralFace = lateralFaceIterator->second;
                   // get crossflow
                   double w = 0.0;//JEH: need to take crossflow from solution vector
                   // compute turbulent crossflow
                   double wt = 0.0;//JEH: need to use turbulent crossflow model
                   // get index of neighboring subchannel
                   std::vector<AMP::Mesh::MeshElement> adjacentCells = d_Mesh->getElementParents(lateralFace, AMP::Mesh::Volume);
                   AMP_INSIST(adjacentCells.size() == 2,"There were not 2 adjacent cells to a lateral gap face");
                   std::vector<double> subchannelCentroid1 = adjacentCells[0].centroid();
                   std::vector<double> subchannelCentroid2 = adjacentCells[1].centroid();
                   size_t subchannelIndex1 = getSubchannelIndex(subchannelCentroid1[0],subchannelCentroid1[1]);
                   size_t subchannelIndex2 = getSubchannelIndex(subchannelCentroid2[0],subchannelCentroid2[1]);
                   size_t neighborSubchannelIndex = 0;
                   AMP::Mesh::MeshElement neighborCell;
                   if (subchannelIndex1 == currentSubchannelIndex) {
                      AMP_INSIST(subchannelIndex2 != currentSubchannelIndex,"Adjacent cells have the same subchannel index.");
                      neighborSubchannelIndex = subchannelIndex2;
                      neighborCell = adjacentCells[1];
                   } else if (subchannelIndex2 == currentSubchannelIndex) {
                      neighborSubchannelIndex = subchannelIndex1;
                      neighborCell = adjacentCells[0];
                   } else {
                      AMP_ERROR("Neither of adjacent cells had the same index as the current subchannel.");
                   }
                   // determine sign of crossflow term
                   double crossflowSign =0.;
                   if (currentSubchannelIndex < neighborSubchannelIndex) {
                      crossflowSign = 1.0;
                   } else if (currentSubchannelIndex > neighborSubchannelIndex) {
                      crossflowSign = -1.0;
                   } else {
                      AMP_ERROR("Adjacent cells have the same subchannel index.");
                   }
                   // get upper and lower axial faces of neighbor cell
                   // ------------------------------------------------
                   std::vector<double> neighborCentroid = neighborCell.centroid();
                   // get axial faces of current cell
                   std::vector<AMP::Mesh::MeshElement> neighborFaces = neighborCell.getElements(AMP::Mesh::Face);
                   AMP::Mesh::MeshElement neighborPlusFace; // upper axial face for current cell
                   AMP::Mesh::MeshElement neighborMinusFace; // lower axial face for current cell
                   // loop over faces of current cell
                   for (std::vector<AMP::Mesh::MeshElement>::iterator face = neighborFaces.begin(); face != neighborFaces.end(); ++face) {
                      std::vector<double> faceCentroid = face->centroid();
                      // if z-coordinates of centroids of the cell and face are not equal,
                      if (!AMP::Utilities::approx_equal(faceCentroid[2],neighborCentroid[2],1.0e-6)) {
                         // if face is above cell centroid,
                         if (faceCentroid[2] > neighborCentroid[2])
                            // face is the upper axial face
                            neighborPlusFace = *face;
                         else
                            // face is the lower axial face
                            neighborMinusFace = *face;
                      }
                   }
      
                   // extract unknowns from solution vector
                   // -------------------------------------
                   std::vector<size_t> neighborPlusDofs;
                   std::vector<size_t> neighborMinusDofs;
                   dof_manager->getDOFs(neighborPlusFace.globalID(),neighborPlusDofs);
                   dof_manager->getDOFs(neighborMinusFace.globalID(),neighborMinusDofs);
                   double m_plus_neighbor = inputVec->getValueByGlobalID(neighborPlusDofs[0]);
                   double h_plus_neighbor = inputVec->getValueByGlobalID(neighborPlusDofs[1]);
                   double p_plus_neighbor = inputVec->getValueByGlobalID(neighborPlusDofs[2]);
                   double m_minus_neighbor = inputVec->getValueByGlobalID(neighborMinusDofs[0]);
                   double h_minus_neighbor = inputVec->getValueByGlobalID(neighborMinusDofs[1]);
                   double p_minus_neighbor = inputVec->getValueByGlobalID(neighborMinusDofs[2]);

                   // compute additional quantities from neighboring cell
                   double m_mid_neighbor = 1.0/2.0*(m_plus_neighbor + m_minus_neighbor);
                   double p_mid_neighbor = 1.0/2.0*(p_plus_neighbor + p_minus_neighbor);

                   // evaluate specific volume at upper face
                   std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_plus_neighbor;
                   volumeArgMap_plus_neighbor.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_plus_neighbor)));
                   volumeArgMap_plus_neighbor.insert(std::make_pair("pressure",new std::vector<double>(1,p_plus_neighbor)));
                   std::vector<double> volumeResult_plus_neighbor(1);
                   d_subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_plus_neighbor,volumeArgMap_plus_neighbor); 
                   double vol_plus_neighbor = volumeResult_plus_neighbor[0];
      
                   // evaluate specific volume at lower face
                   std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_minus_neighbor;
                   volumeArgMap_minus_neighbor.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_minus_neighbor)));
                   volumeArgMap_minus_neighbor.insert(std::make_pair("pressure",new std::vector<double>(1,p_minus_neighbor)));
                   std::vector<double> volumeResult_minus_neighbor(1);
                   d_subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_minus_neighbor,volumeArgMap_minus_neighbor); 
                   double vol_minus_neighbor = volumeResult_minus_neighbor[0];

                   double h_axialDonor_neighbor;
                   double vol_axialDonor_neighbor;
                   if (m_mid_neighbor >= 0.0) {
                      h_axialDonor_neighbor = h_minus_neighbor;
                      vol_axialDonor_neighbor = vol_minus_neighbor;
                   } else {
                      h_axialDonor_neighbor = h_plus_neighbor;
                      vol_axialDonor_neighbor = vol_plus_neighbor;
                   }

                   double rho_mid_neighbor = 1.0/vol_axialDonor_neighbor;
                   double u_mid_neighbor = m_mid_neighbor*vol_axialDonor_neighbor/area;
                      
                   double h_lateralDonor;
                   double u_lateralDonor;
                   if (crossflowSign*w >= 0.0) {
                      h_lateralDonor = h_axialDonor;
                      u_lateralDonor = u_mid;
                   } else {
                      h_lateralDonor = h_axialDonor_neighbor;
                      u_lateralDonor = u_mid_neighbor;
                   }
             
                   // evaluate temperature for neighbor cell
                   std::map<std::string, boost::shared_ptr<std::vector<double> > > temperatureArgMap_mid_neighbor;
                   temperatureArgMap_mid_neighbor.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_axialDonor_neighbor)));
                   temperatureArgMap_mid_neighbor.insert(std::make_pair("pressure",new std::vector<double>(1,p_mid_neighbor)));
                   std::vector<double> temperatureResult_mid_neighbor(1);
                   d_subchannelPhysicsModel->getProperty("Temperature",temperatureResult_mid_neighbor,temperatureArgMap_mid_neighbor); 
                   double T_mid_neighbor = temperatureResult_mid_neighbor[0];
      
                   // evaluate conductivity for cell
                   std::map<std::string, boost::shared_ptr<std::vector<double> > > conductivityArgMap_mid_neighbor;
                   conductivityArgMap_mid_neighbor.insert(std::make_pair("conductivity",new std::vector<double>(1,T_mid_neighbor)));
                   conductivityArgMap_mid_neighbor.insert(std::make_pair("density",new std::vector<double>(1,rho_mid_neighbor)));
                   std::vector<double> conductivityResult_mid_neighbor(1);
                   d_subchannelPhysicsModel->getProperty("ThermalConductivity",conductivityResult_mid_neighbor,conductivityArgMap_mid_neighbor); 
                   double k_mid_neighbor = conductivityResult_mid_neighbor[0];

                   // compute thermal conductivity across gap
                   double k_gap = 2.0*k_mid*k_mid_neighbor/(k_mid + k_mid_neighbor);

                   // compute distance between centroids of cells adjacent to gap
                   double lx = std::abs(neighborCentroid[0] - cellCentroid[0]);
                   double ly = std::abs(neighborCentroid[0] - cellCentroid[0]);
                   double l = std::max(lx,ly);

                   // compute gap width
                   double s = d_pitch - d_diameter;//JEH: need to get information from mesh

                   double conductance = 1.0*k_gap/l;

                   // add to sums
                   mass_crossflow_sum += crossflowSign*w;
                   energy_crossflow_sum += crossflowSign*w*h_lateralDonor;
                   energy_turbulence_sum += wt*(h_axialDonor - h_axialDonor_neighbor);
                   energy_conduction_sum += conductance*s*(T_mid - T_mid_neighbor);
                   axial_crossflow_sum += crossflowSign*w*u_lateralDonor;
                      
                }// end if (lateralFaceIterator != lateralFaceMap.end()) {
             }// end loop over gap faces

             // loop over rods

             // calculate residuals for current cell
             // ------------------------------------
             double R_m = m_plus - m_minus
                        + mass_crossflow_sum; // mass
             double R_h = (m_plus*h_plus - m_minus*h_minus)/dz
                        + energy_crossflow_sum/dz
                        - energy_heatflux_sum
                        + energy_turbulence_sum/dz
                        + energy_conduction_sum
                        - energy_direct_heating_sum; // energy
             double R_p = m_plus*u_plus - m_minus*u_minus
                        + axial_crossflow_sum
                        + area*(p_plus-p_minus)
                        + g*area*dz*std::cos(d_theta)/vol_axialDonor
                        + 1.0/(2.0*area)*(dz*d_friction/D + d_K)*std::abs(m_mid)*m_mid*vol_axialDonor
                        + axial_turbulence_sum; // axial momentum

             // put residuals into global residual vector
             outputVec->setValueByGlobalID(plusDofs[0], R_m);
             outputVec->setValueByGlobalID(plusDofs[1], R_h);
             outputVec->setValueByGlobalID(minusDofs[2], R_p);

             // impose boundary conditions
             // --------------------------
             // if face is exit face,
             if (AMP::Utilities::approx_equal(plusFaceCentroid[2],height)){
                // impose fixed exit pressure boundary condition
                outputVec->setValueByGlobalID(plusDofs[2],p_plus-d_Pout);
             }
             if (AMP::Utilities::approx_equal(minusFaceCentroid[2],0.0)){
                // impose fixed inlet axial mass flow rate boundary condition
                outputVec->setValueByGlobalID(minusDofs[0],m_minus-d_min);

                // evaluate enthalpy at inlet
                std::map<std::string, boost::shared_ptr<std::vector<double> > > enthalpyArgMap_inlet;
                enthalpyArgMap_inlet.insert(std::make_pair("temperature",new std::vector<double>(1,d_Tin)));
                enthalpyArgMap_inlet.insert(std::make_pair("pressure",new std::vector<double>(1,p_minus)));
                std::vector<double> enthalpyResult_inlet(1);
                d_subchannelPhysicsModel->getProperty("Enthalpy",enthalpyResult_inlet,enthalpyArgMap_inlet); 
                double h_eval = enthalpyResult_inlet[0];
                // impose fixed inlet temperature boundary condition
                outputVec->setValueByGlobalID(minusDofs[0],h_minus-h_eval);
             }
          }// end loop over cells of current subchannel
        }// end if(d_ownSubchannel[isub])
      }// end loop over subchannels

      // loop over lateral faces
      /*
      AMP::Mesh::MeshIterator face = d_Mesh->getIterator(AMP::Mesh::Face, 0); // iterator for cells of mesh
      for (; face != face.end(); face++) {
         std::vector<double> faceCentroid = face->centroid();
         std::map<std::vector<double>,AMP::Mesh::MeshElement>::iterator lateralFaceIterator = lateralFaceMap.find(faceCentroid);
         if (lateralFaceIterator != lateralFaceMap.end()) {
            // get face
            AMP::Mesh::MeshElement lateralFace = lateralFaceIterator->second;
            // get crossflow from solution vector
            std::vector<size_t> gapDofs;
            dof_manager->getDOFs(lateralFace.globalID(),gapDofs);
            double w = inputVec->getValueByGlobalID(gapDofs[0]);

            // get adjacent cells
            std::vector<AMP::Mesh::MeshElement> adjacentCells = d_Mesh->getElementParents(lateralFace, AMP::Mesh::Volume);
            AMP_INSIST(adjacentCells.size() == 2,"There were not 2 adjacent cells to a lateral gap face");

            double R_w = 0.0;
            outputVec->setValueByGlobalID(gapDofs[0],R_w);
         }
      }// end loop over lateral faces
*/

}// end of apply function

boost::shared_ptr<OperatorParameters> SubchannelFourEqNonlinearOperator :: 
getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u) 
{
  boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));

  tmp_db->putString("name","SubchannelFourEqNonlinearOperator");

  boost::shared_ptr<SubchannelOperatorParameters> outParams(new SubchannelOperatorParameters(tmp_db));
  return outParams;
}

AMP::LinearAlgebra::Vector::shared_ptr  SubchannelFourEqNonlinearOperator::subsetInputVector(AMP::LinearAlgebra::Vector::shared_ptr vec)
{
  AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
  // Subset the vectors, they are simple vectors and we need to subset for the current comm instead of the mesh
  if(d_Mesh.get() != NULL) {
    AMP::LinearAlgebra::VS_Comm commSelector( var->getName(), d_Mesh->getComm() );
    AMP::LinearAlgebra::Vector::shared_ptr commVec = vec->select(commSelector, var->getName());
    return commVec->subsetVectorForVariable(var);
  } else {
    return vec->subsetVectorForVariable(var);
  }
}

AMP::LinearAlgebra::Vector::const_shared_ptr  SubchannelFourEqNonlinearOperator::subsetInputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec)
{
  AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
  // Subset the vectors, they are simple vectors and we need to subset for the current comm instead of the mesh
  if(d_Mesh.get() != NULL) {
    AMP::LinearAlgebra::VS_Comm commSelector( var->getName(), d_Mesh->getComm() );
    AMP::LinearAlgebra::Vector::const_shared_ptr commVec = vec->constSelect(commSelector, var->getName());
    return commVec->constSubsetVectorForVariable(var);
  } else {
    return vec->constSubsetVectorForVariable(var);
  }
}

AMP::LinearAlgebra::Vector::shared_ptr  SubchannelFourEqNonlinearOperator::subsetOutputVector(AMP::LinearAlgebra::Vector::shared_ptr vec)
{
  AMP::LinearAlgebra::Variable::shared_ptr var = getOutputVariable();
  // Subset the vectors, they are simple vectors and we need to subset for the current comm instead of the mesh
  if(d_Mesh.get() != NULL) {
    AMP::LinearAlgebra::VS_Comm commSelector( var->getName(), d_Mesh->getComm() );
    AMP::LinearAlgebra::Vector::shared_ptr commVec = vec->select(commSelector, var->getName());
    return commVec->subsetVectorForVariable(var);
  } else {
    return vec->subsetVectorForVariable(var);
  }
}

AMP::LinearAlgebra::Vector::const_shared_ptr  SubchannelFourEqNonlinearOperator::subsetOutputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec)
{
  AMP::LinearAlgebra::Variable::shared_ptr var = getOutputVariable();
    // Subset the vectors, they are simple vectors and we need to subset for the current comm instead of the mesh
    if(d_Mesh.get() != NULL) {
        AMP::LinearAlgebra::VS_Comm commSelector( var->getName(), d_Mesh->getComm() );
        AMP::LinearAlgebra::Vector::const_shared_ptr commVec = vec->constSelect(commSelector, var->getName());
        return commVec->constSubsetVectorForVariable(var);
    } else {
        return vec->constSubsetVectorForVariable(var);
    }
}

}
}
