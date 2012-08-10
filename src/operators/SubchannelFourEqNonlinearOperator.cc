
#include "operators/SubchannelFourEqNonlinearOperator.h"
#include "operators/SubchannelOperatorParameters.h"
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
      d_min  = getDoubleParameter(myparams,"Axial_Flow_Rate",0.3522);  
      d_win  = getDoubleParameter(myparams,"Cross_Flow_Rate",0.3522);  
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

// apply
void SubchannelFourEqNonlinearOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
    AMP::LinearAlgebra::Vector::shared_ptr &r, const double a, const double b)
{

      // ensure that solution and residual vectors aren't NULL
      AMP_INSIST( ((r.get()) != NULL), "NULL Residual Vector" );
      AMP_INSIST( ((u.get()) != NULL), "NULL Solution Vector" );

      // calculate extra parameters
      const double pi = 4.0*atan(1.0); // pi
      const double g = 9.805;          // acceleration due to gravity [m/s2]
      // assuming square pitch
      const double perimeter = 4.0*(d_pitch-d_diameter) + pi*d_diameter;    // wetted perimeter
      const double A = std::pow(d_pitch,2) - pi*std::pow(d_diameter,2)/4.0; // flow area
      const double D = 4.0*A/perimeter;                                     // hydraulic diameter

      // Subset the vectors
      AMP::LinearAlgebra::Vector::shared_ptr inputVec = subsetInputVector( u );
      AMP::LinearAlgebra::Vector::shared_ptr outputVec = subsetOutputVector( r );

      AMP::Discretization::DOFManager::shared_ptr dof_manager = inputVec->getDOFManager();

      // Iterate over the subchannel and get the meshes
      int subChannel = 0;

      AMP::Mesh::MeshIterator cell = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
      AMP::Mesh::MeshIterator end_cell = cell.end();

      // get interval lengths from mesh
      std::vector<double> box = d_Mesh->getBoundingBox();
      const double min_z = box[4];
      const double max_z = box[5];
      const double height = max_z - min_z;

      ///assuming singl subchannel
      const int numCells = cell.size();
      const int numFaces = numCells + 1;

      // compute element heights
      std::vector<double> del_z(numCells);
      // assuming uniform mesh
      for (int j=0; j<numCells; j++) {
        del_z[j] = height/numCells; 
      }

      // create vector of axial positions
      // axial positions are only used if some rod power shape is assumed
      std::vector<double> z(numFaces);
      z[0] = 0.0;
      for( int j=1; j<numFaces; j++) {
        z[j] = z[j-1] + del_z[j-1];
      } 
      
      // calculate residual for axial momentum equations
      double R_h, R_pa, R_pl, R_m; 
      int j = 1;
      std::vector<double> zcoordIt, xcoordIt;
      for(size_t icell = 0; icell < cell.size(); ++icell, ++j, ++cell){

        std::vector<AMP::Mesh::MeshElement> d_currFaces = cell->getElements(AMP::Mesh::Face);

        std::vector<double> cellCentroid = cell->centroid();

        std::map<double,AMP::Mesh::MeshElement> xyFace, gapFace;
        for(size_t idx=0; idx<d_currFaces.size(); idx++ )
        {
          std::vector<double> ctrd = d_currFaces[idx].centroid();
          std::vector<AMP::Mesh::MeshElement> nodes = d_currFaces[idx].getElements(AMP::Mesh::Vertex);

          bool is_valid = true;
          for (size_t j=0; j<nodes.size(); ++j) {
            std::vector<double> coord = nodes[j].coord();
            if ( !AMP::Utilities::approx_equal(coord[2],ctrd[2], 1e-6) )
              is_valid = false;
          }
          if ( is_valid ) {
            xyFace.insert(std::pair<double,AMP::Mesh::MeshElement>(ctrd[2],d_currFaces[idx]));
            zcoordIt.push_back(ctrd[2]);
          }else if((d_Mesh->getElementParents(d_currFaces[idx], AMP::Mesh::Volume)).size() > 1 ){
              gapFace.insert(std::pair<double,AMP::Mesh::MeshElement>(ctrd[0],d_currFaces[idx]));
              xcoordIt.push_back(ctrd[2]);
          }
        }

        std::vector<size_t> axialDofs;
        std::vector<double> hAxial(xyFace.size()), PAxial(xyFace.size()), mAxial(xyFace.size());
        std::vector<double> wGap(gapFace.size());

        for(size_t idxy=0; idxy<xyFace.size(); ++idxy){
          dof_manager->getDOFs( (xyFace.find(zcoordIt[idxy])->second).globalID(), axialDofs);

          hAxial[idxy]  = inputVec->getValueByGlobalID(axialDofs[0]);
          PAxial[idxy]  = inputVec->getValueByGlobalID(axialDofs[1]); 
          mAxial[idxy]  = inputVec->getValueByGlobalID(axialDofs[2]); 
        }

        std::vector<size_t> gapDofs;
        for(size_t idxy=0; idxy<gapFace.size(); ++idxy){
          dof_manager->getDOFs( (gapFace.find(xcoordIt[idxy])->second).globalID(), gapDofs);
          wGap[idxy]   = inputVec->getValueByGlobalID(gapDofs[0]);
        }

        std::vector<bool> isNeighbor(gapFace.size(), false);
        for(size_t idxy=0; idxy<gapFace.size(); ++idxy){
          std::vector<AMP::Mesh::MeshElement> faceParents = (d_Mesh->getElementParents((gapFace.find(xcoordIt[idxy])->second), AMP::Mesh::Volume));
          std::vector<double> currCentroid = faceParents[0].centroid();
            if ( !AMP::Utilities::approx_equal(cellCentroid[0], currCentroid[0], 1e-6) && !AMP::Utilities::approx_equal(cellCentroid[1], currCentroid[1], 1e-6)){
               isNeighbor[0] = true;
            }else{
               isNeighbor[1] = true;
            }
        }

        //lateral momentuem equation
        for(size_t idxy=0; idxy<gapFace.size(); ++idxy){
          R_pl = wGap[idxy]- d_win;
          outputVec->setValueByGlobalID(gapDofs[idxy], R_pl);
        }

        //axial momentum equation
        double h_avg   = (1.0/2.0)*(hAxial[0] + hAxial[1]); // enthalpy evaluated at cell center
        double p_avg   = (1.0/2.0)*(PAxial[0] + PAxial[1]);       // pressure evaluated at cell center
        double m_avg   = (1.0/2.0)*(mAxial[0] + mAxial[1]);       // pressure evaluated at cell center

        // evaluate density at upper face
        std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_plus;
        volumeArgMap_plus.insert(std::make_pair("enthalpy",new std::vector<double>(1,hAxial[1])));
        volumeArgMap_plus.insert(std::make_pair("pressure",new std::vector<double>(1,PAxial[1])));
        std::vector<double> volumeResult_plus(1);
        d_subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_plus,volumeArgMap_plus); 
        double rho_plus = 1.0/volumeResult_plus[0];

        // evaluate density at lower face
        std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_minus;
        volumeArgMap_minus.insert(std::make_pair("enthalpy",new std::vector<double>(1,hAxial[0])));
        volumeArgMap_minus.insert(std::make_pair("pressure",new std::vector<double>(1,PAxial[0])));
        std::vector<double> volumeResult_minus(1);
        d_subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_minus,volumeArgMap_minus); 
        double rho_minus = 1.0/volumeResult_minus[0];

        // evaluate density at cell center
        std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_avg;
        volumeArgMap_avg.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_avg)));
        volumeArgMap_avg.insert(std::make_pair("pressure",new std::vector<double>(1,p_avg)));
        std::vector<double> volumeResult_avg(1);
        d_subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_avg,volumeArgMap_avg);
        double rho_avg = 1.0/volumeResult_avg[0];

        double u_plus  = mAxial[1]/ (A*rho_plus);  // velocity evaluated at upper face
        double u_minus = mAxial[0]/ (A*rho_minus); // velocity evaluated at lower face

        // evaluate residual: axial momentum equation
        R_pa = (mAxial[1]*u_plus - mAxial[0]*u_minus)
          + g * d_area[subChannel] * del_z[j-1] * rho_avg * std::cos(d_theta)
          + (1.0/(2.0*d_area[subChannel]))*(del_z[j-1] * (d_friction/D) + d_K)* std::abs(m_avg)*(m_avg/rho_avg)
          + d_area[subChannel] * (PAxial[1]- PAxial[0]);



      }//end for cell

}

/*
// apply
void SubchannelFourEqNonlinearOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
AMP::LinearAlgebra::Vector::shared_ptr &r, const double a, const double b)
{

// ensure that solution and residual vectors aren't NULL
AMP_INSIST( ((r.get()) != NULL), "NULL Residual Vector" );
AMP_INSIST( ((u.get()) != NULL), "NULL Solution Vector" );

// calculate extra parameters
const double pi = 4.0*atan(1.0); // pi
const double g = 9.805;          // acceleration due to gravity [m/s2]
// assuming square pitch
const double perimeter = 4.0*(d_pitch-d_diameter) + pi*d_diameter;    // wetted perimeter
const double A = std::pow(d_pitch,2) - pi*std::pow(d_diameter,2)/4.0; // flow area
const double D = 4.0*A/perimeter;                                     // hydraulic diameter

// Subset the vectors
AMP::LinearAlgebra::Vector::shared_ptr inputVec = subsetInputVector( u );
AMP::LinearAlgebra::Vector::shared_ptr outputVec = subsetOutputVector( r );

AMP::Discretization::DOFManager::shared_ptr dof_manager = inputVec->getDOFManager();

// get the Iterators for the subchannel mesh
AMP::Mesh::MeshIterator begin_face = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(d_Mesh, 0);
AMP::Mesh::MeshIterator end_face   = begin_face.end();


get boundary values: u has ordering:
\f[ \vec{u}=\left[\begin{array}{c}
h_{0^+}\\
p_{0^+}\\
\vdots\\
h_{J^+}\\
p_{J^+}\\
\end{array}\right] \f]

const int numFaces = begin_face.size() ;
const int numCells = numFaces - 1;

std::vector<size_t> dofs;
dof_manager->getDOFs( begin_face->globalID(), dofs );

double h_in  = inputVec->getValueByGlobalID(dofs[0]);
double P_in  = inputVec->getValueByGlobalID(dofs[1]);    

// evaluate enthalpy at inlet
std::map<std::string, boost::shared_ptr<std::vector<double> > > enthalpyArgMap;
enthalpyArgMap.insert(std::make_pair("temperature",new std::vector<double>(1,d_Tin)));
enthalpyArgMap.insert(std::make_pair("pressure",   new std::vector<double>(1,P_in)));
std::vector<double> enthalpyResult(1);
d_subchannelPhysicsModel->getProperty("Enthalpy",enthalpyResult,enthalpyArgMap); 
double h_eval = enthalpyResult[0];

// get interval lengths from mesh
std::vector<double> box = d_Mesh->getBoundingBox();
const double min_z = box[4];
const double max_z = box[5];
const double height = max_z - min_z;

// compute element heights
std::vector<double> del_z(numCells);
// assuming uniform mesh
for (int j=0; j<numCells; j++) {
del_z[j] = height/numCells; 
}

// create vector of axial positions
// axial positions are only used if some rod power shape is assumed
std::vector<double> z(numFaces);
z[0] = 0.0;
for( int j=1; j<numFaces; j++) {
  z[j] = z[j-1] + del_z[j-1];
} 

// compute the enthalpy change in each interval
std::vector<double> dh(numCells);
if (d_source == "averageCladdingTemperature") {
  AMP_ERROR("Heat source type 'averageCladdingTemperature' not yet implemented.");
} else if (d_source == "averageHeatFlux") {
  AMP_ERROR("Heat source type 'averageHeatFlux' not yet implemented.");
} else if (d_source == "totalHeatGeneration") {
  if (d_heatShape == "Sinusoidal") {
    // sinusoidal
    for (int j=0; j<numCells; j++){
      double flux = d_Q/(2.0*pi*d_diameter*del_z[j]) * (std::cos(pi*z[j]/height) - std::cos(pi*z[j+1]/height));
      double lin = d_Q/(2.0*del_z[j])                * (std::cos(pi*z[j]/height) - std::cos(pi*z[j+1]/height));
      double flux_sum = 4.0*pi*d_diameter*1.0/4.0*flux;
      double lin_sum = 4.0*d_gamma*1.0/4.0*lin;
      dh[j] = del_z[j] / d_m * (flux_sum + lin_sum);
    }
  } else {
    AMP_ERROR("Heat shape '"+d_heatShape+" is invalid");
  }
} else {
  AMP_ERROR("Heat source type '"+d_source+"' is invalid");
}

// calculate residual for axial momentum equations
double R_h, R_p; 
int j = 1;
AMP::Mesh::MeshIterator face = begin_face;
for(size_t iface = 0; iface < begin_face.size(); ++iface, ++j){
  // ======================================================
  // energy residual
  // ======================================================
  if (face == begin_face){

    evaluate first residual entry, corresponding to inlet enthalpy:
      \f[ R_0 = h_{in} - h(T_{in},p_{1-})\f]

      R_h = h_in - h_eval;
  } else {
    // residual at face corresponds to cell below
    dof_manager->getDOFs( face->globalID(), dofs );
    double h_plus   = inputVec->getValueByGlobalID(dofs[0]); // enthalpy evaluated at lower face
    --face;
    dof_manager->getDOFs( face->globalID(), dofs );
    double h_minus  = inputVec->getValueByGlobalID(dofs[0]); // enthalpy evaluated at lower face
    ++face;

    R_h = h_plus - h_minus - dh[j-2];
  }

  // ======================================================
  // axial momentum residual
  // ======================================================
  // residual at face corresponds to cell above
  dof_manager->getDOFs( face->globalID(), dofs );
  double h_minus = inputVec->getValueByGlobalID(dofs[0]); // enthalpy evaluated at lower face
  double p_minus = inputVec->getValueByGlobalID(dofs[1]); // pressure evaluated at lower face
  if (face == end_face - 1){
    R_p = p_minus - d_Pout;
  } else {
    ++face;
    dof_manager->getDOFs( face->globalID(), dofs );
    double h_plus  = inputVec->getValueByGlobalID(dofs[0]); // enthalpy evaluated at upper face
    double p_plus  = inputVec->getValueByGlobalID(dofs[1]); // pressure evaluated at upper face
    --face;

    double h_avg   = (1.0/2.0)*(h_minus + h_plus); // enthalpy evaluated at cell center
    double p_avg   = (1.0/2.0)*(p_minus + p_plus);       // pressure evaluated at cell center

    // evaluate density at upper face
    std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_plus;
    volumeArgMap_plus.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_plus)));
    volumeArgMap_plus.insert(std::make_pair("pressure",new std::vector<double>(1,p_plus)));
    std::vector<double> volumeResult_plus(1);
    d_subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_plus,volumeArgMap_plus); 
    double rho_plus = 1.0/volumeResult_plus[0];

    // evaluate density at lower face
    std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_minus;
    volumeArgMap_minus.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_minus)));
    volumeArgMap_minus.insert(std::make_pair("pressure",new std::vector<double>(1,p_minus)));
    std::vector<double> volumeResult_minus(1);
    d_subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_minus,volumeArgMap_minus); 
    double rho_minus = 1.0/volumeResult_minus[0];

    // evaluate density at cell center
    std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_avg;
    volumeArgMap_avg.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_avg)));
    volumeArgMap_avg.insert(std::make_pair("pressure",new std::vector<double>(1,p_avg)));
    std::vector<double> volumeResult_avg(1);
    d_subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_avg,volumeArgMap_avg);
    double rho_avg = 1.0/volumeResult_avg[0];

    double u_plus  = d_m / (A*rho_plus);  // velocity evaluated at upper face
    double u_minus = d_m / (A*rho_minus); // velocity evaluated at lower face

    // evaluate residual: axial momentum equation
    R_p = (d_m/A)*(u_plus - u_minus)
      + g * del_z[j-1] * rho_avg * std::cos(d_theta)
      + (1.0/2.0)*(del_z[j-1] * d_friction/D + d_K)* std::abs(d_m/(A*rho_avg))*(d_m/A)
      + p_plus - p_minus;
  }

  // put residual value in residual vector
  dof_manager->getDOFs( face->globalID(), dofs );
  outputVec->setValueByGlobalID(dofs[0], R_h);
  outputVec->setValueByGlobalID(dofs[1], R_p);
  ++face;
}

if(f.get() == NULL) {
  outputVec->scale(a);
} else {
  AMP::LinearAlgebra::Vector::shared_ptr fInternal = subsetInputVector( f );
  if(fInternal.get() == NULL) {
    outputVec->scale(a);
  } else {
    outputVec->axpby(b, a, fInternal);
  }
}
}
*/

// JEH: what is the purpose of this function?
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
