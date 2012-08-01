
#include "operators/SubchannelTwoEqNonlinearOperator.h"
#include "operators/SubchannelOperatorParameters.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

#include <string>


namespace AMP {
namespace Operator {

// reset
void SubchannelTwoEqNonlinearOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
{
      boost::shared_ptr<SubchannelOperatorParameters> myparams = 
        boost::dynamic_pointer_cast<SubchannelOperatorParameters>(params);

      AMP_INSIST( ((myparams.get()) != NULL), "NULL parameters" );
      AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );

      d_Pout = getDoubleParameter(myparams,"Exit_Pressure",15.5132e6);
      d_Tin  = getDoubleParameter(myparams,"Inlet_Temperature",569.26);  
      d_m    = getDoubleParameter(myparams,"Mass_Flow_Rate",0.3522);  
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
      }

      d_subchannelPhysicsModel = myparams->d_subchannelPhysicsModel;  
}

// function used in reset to get double parameter or set default if missing
double SubchannelTwoEqNonlinearOperator::getDoubleParameter(	boost::shared_ptr<SubchannelOperatorParameters> myparams,
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
std::string SubchannelTwoEqNonlinearOperator::getStringParameter(	boost::shared_ptr<SubchannelOperatorParameters> myparams,
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
void SubchannelTwoEqNonlinearOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
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
      AMP::Mesh::MeshIterator face     = d_Mesh->getIterator(AMP::Mesh::Face, 0);
      AMP::Mesh::MeshIterator end_face = face.end();

      /**
        get boundary values: u has ordering:
        \f[ \vec{u}=\left[\begin{array}{c}
                          h_{in}\\
                          p_{1^-}\\
                          \vdots\\
                          p_{J^+}\\
                          \end{array}\right] \f]
        */
      const int numCells = face.size() - 1;
      const int numFaces = face.size() ;

      std::vector<size_t> dofs;
      dof_manager->getDOFs( face->globalID(), dofs );

      double h_in  = inputVec->getValueByGlobalID(dofs[0]);
      double P_in  = inputVec->getValueByGlobalID(dofs[1]);    

      // evaluate enthalpy at inlet
      std::map<std::string, boost::shared_ptr<std::vector<double> > > enthalpyArgMap;
      enthalpyArgMap.insert(std::make_pair("temperature",new std::vector<double>(1,d_Tin)));
      enthalpyArgMap.insert(std::make_pair("pressure",   new std::vector<double>(1,P_in)));
      std::vector<double> enthalpyResult(1);
      d_subchannelPhysicsModel->getProperty("Enthalpy",enthalpyResult,enthalpyArgMap); 
      double h_eval = enthalpyResult[0];

      /**
        evaluate first residual entry, corresponding to inlet enthalpy:
        \f[ R_0 = h_{in} - h(T_{in},p_{1-})\f]
        */
      double R_b = h_in - h_eval;

      // put first residual value into residual vector
      outputVec->setValueByGlobalID(dofs[0], R_b);

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
          // assuming cosine power shape
          for (int j=0; j<numCells; j++){
              double flux = d_Q/(2.0*pi*d_diameter*del_z[j]) * (std::cos(pi*z[j]/height) - std::cos(pi*z[j+1]/height));
              double lin = d_Q/(2.0*del_z[j])                * (std::cos(pi*z[j]/height) - std::cos(pi*z[j+1]/height));
              double flux_sum = 4.0*pi*d_diameter*1.0/4.0*flux;
              double lin_sum = 4.0*d_gamma*1.0/4.0*lin;
              dh[j] = del_z[j] / d_m * (flux_sum + lin_sum);
          }
      } else {
          AMP_ERROR("Heat source type '"+d_source+"' is invalid");
      }

      // strongly impose outlet pressure boundary condition
      dof_manager->getDOFs( end_face->globalID(), dofs );
      inputVec->setValueByGlobalID(dofs[1], d_Pout);
      outputVec->setValueByGlobalID(dofs[1], 0.0);

 
      // calculate residual for axial momentum equations
      double h_minus = h_in;
      double h_plus = h_in;
      int j = 1;
      for( ; face != end_face; ++j){
          h_minus = h_plus;                            // enthalpy evaluated at lower face
          h_plus  = h_minus + dh[j-1];                 // enthalpy evaluated at upper face
          double h_avg   = 1.0/2.0*(h_minus + h_plus); // enthalpy evaluated at cell center

          dof_manager->getDOFs( face->globalID(), dofs );
          double p_minus = inputVec->getValueByGlobalID(dofs[1]);   // pressure evaluated at lower face

          ++face;
          dof_manager->getDOFs( face->globalID(), dofs );
          double p_plus  = inputVec->getValueByGlobalID(dofs[1]); // pressure evaluated at upper face

          double p_avg   = 1.0/2.0*(p_minus + p_plus);       // pressure evaluated at cell center

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
          R_b = (d_m/A)*(u_plus - u_minus)
              + g * del_z[j-1] * rho_avg * std::cos(d_theta) + 
              + 1.0/2.0*(del_z[j-1] * d_friction/D + d_K)* std::abs(d_m/(A*rho_avg))*(d_m/A)
              + p_plus - p_minus;

          // put residual value in residual vector
          outputVec->setValueByGlobalID(dofs[1], R_b);
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

// JEH: what is the purpose of this function?
boost::shared_ptr<OperatorParameters> SubchannelTwoEqNonlinearOperator :: 
    getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u) 
{
        boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));

        tmp_db->putString("name","SubchannelTwoEqNonlinearOperator");

        boost::shared_ptr<SubchannelOperatorParameters> outParams(new SubchannelOperatorParameters(tmp_db));
        return outParams;
}

AMP::LinearAlgebra::Vector::shared_ptr  SubchannelTwoEqNonlinearOperator::subsetInputVector(AMP::LinearAlgebra::Vector::shared_ptr vec)
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

AMP::LinearAlgebra::Vector::const_shared_ptr  SubchannelTwoEqNonlinearOperator::subsetInputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec)
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

AMP::LinearAlgebra::Vector::shared_ptr  SubchannelTwoEqNonlinearOperator::subsetOutputVector(AMP::LinearAlgebra::Vector::shared_ptr vec)
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

AMP::LinearAlgebra::Vector::const_shared_ptr  SubchannelTwoEqNonlinearOperator::subsetOutputVector(AMP::LinearAlgebra::Vector::const_shared_ptr vec)
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
