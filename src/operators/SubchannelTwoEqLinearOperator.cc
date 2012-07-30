
#include "operators/SubchannelTwoEqLinearOperator.h"
#include "operators/SubchannelOperatorParameters.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

#include <string>


namespace AMP {
namespace Operator {

// reset
void SubchannelTwoEqLinearOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
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

      // get subchannel physics model
      d_subchannelPhysicsModel = myparams->d_subchannelPhysicsModel;  

      // get frozen solution
      if ((myparams->d_frozenSolution.get()) != NULL) d_frozenVec = myparams->d_frozenSolution;

      // get option of printing Jacobian
      if ((myparams->d_db)->keyExists("Print_Jacobian")) d_printJacobian = (myparams->d_db)->getBoolWithDefault("Print_Jacobian",false);
}

// function used in reset to get double parameter or set default if missing
double SubchannelTwoEqLinearOperator::getDoubleParameter(	boost::shared_ptr<SubchannelOperatorParameters> myparams,
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
std::string SubchannelTwoEqLinearOperator::getStringParameter(	boost::shared_ptr<SubchannelOperatorParameters> myparams,
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
void SubchannelTwoEqLinearOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
    AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a, const double b)
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

      AMP_INSIST( (d_frozenVec.get() != NULL), "Null Frozen Vector inside Jacobian" );

      // get solution sizes
      const size_t numCells = inputVec->getLocalSize() - 2;
      const size_t numFaces = numCells + 1;
      d_solutionSize = numFaces + 1;

      // create Jacobian matrix if requested
      if (d_printJacobian){
         // allocate memory for matrix
         d_Jacobian = new double*[numFaces+1];
         for (size_t i = 0; i <= numFaces; i++)
            d_Jacobian[i] = new double[numFaces+1];
         // initialize entries to zero
         for (size_t i = 0; i <= numFaces; i++){
            for (size_t j = 0; j <= numFaces; j++)
               d_Jacobian[i][j] = 0.0;
         }
      }

      // get interval lengths from mesh
      std::vector<double> box = d_Mesh->getBoundingBox();
      const double min_z = box[4];
      const double max_z = box[5];
      const double height = max_z - min_z;

      // compute element heights
      std::vector<double> del_z(numCells);
      // assuming uniform mesh
      for (size_t j=0; j<numCells; j++) {
        del_z[j] = height/numCells; 
      }

      // create vector of axial positions
      // axial positions are only used if some rod power shape is assumed
      std::vector<double> z(numFaces);
      z[0] = 0.0;
      for( size_t j=1; j<numFaces; j++) {
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
          for (size_t j=0; j<numCells; j++){
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
      outputVec->setValueByLocalID(numFaces, 0.0);

      // compute Jacobian entry J[0][0]
      double J_0_0 = 1.0;
      if (d_printJacobian) d_Jacobian[0][0] = J_0_0;

      // compute Jacobian entry J[0][1]
      double p_in = d_frozenVec->getValueByLocalID(1);
      double J_0_1 = -1.0*dhdp(d_Tin,p_in);
      if (d_printJacobian) d_Jacobian[0][1] = J_0_1;

      // get delta solution variables
      double du_0 = inputVec->getValueByLocalID(0); // du[0]
      double du_1 = inputVec->getValueByLocalID(1); // du[1]

      // compute entry 0 of Jacobian matrix-vector product J*du
      double R_0 = J_0_0*du_0 + J_0_1*du_1;

      // put residual value in residual vector
      outputVec->setValueByLocalID(0, R_0);

      // calculate residual for axial momentum equations
      double h_in = d_frozenVec->getValueByLocalID(0);
      double h_minus = h_in;
      double h_plus = h_in;
      for( size_t j=1; j <= numFaces-1; j++) {
          h_minus = h_plus;             // enthalpy evaluated at lower face
          h_plus  = h_minus + dh[j-1];  // enthalpy evaluated at upper face

          double p_plus  = d_frozenVec->getValueByLocalID(j+1); // pressure evaluated at upper face
          double p_minus = d_frozenVec->getValueByLocalID(j);   // pressure evaluated at lower face

          // evaluate specific volume at upper face
          std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_plus;
          volumeArgMap_plus.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_plus)));
          volumeArgMap_plus.insert(std::make_pair("pressure",new std::vector<double>(1,p_plus)));
          std::vector<double> volumeResult_plus(1);
          d_subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_plus,volumeArgMap_plus); 
          double v_plus = volumeResult_plus[0];

          // evaluate specific volume at lower face
          std::map<std::string, boost::shared_ptr<std::vector<double> > > volumeArgMap_minus;
          volumeArgMap_minus.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_minus)));
          volumeArgMap_minus.insert(std::make_pair("pressure",new std::vector<double>(1,p_minus)));
          std::vector<double> volumeResult_minus(1);
          d_subchannelPhysicsModel->getProperty("SpecificVolume",volumeResult_minus,volumeArgMap_minus); 
          double v_minus = volumeResult_minus[0];

          // evaluate derivatives
          double dvdh_plus  = dvdh(h_plus, p_plus);
          double dvdh_minus = dvdh(h_minus,p_minus);
          double dvdp_plus  = dvdp(h_plus, p_plus);
          double dvdp_minus = dvdp(h_minus,p_minus);
      
          // compute Jacobian entry J[j][0]
          double J_j_0 = std::pow(d_m/A,2)*(dvdh_plus - dvdh_minus) - 2.0*g*del_z[j-1]*std::cos(d_theta)*
             (dvdh_plus + dvdh_minus)/std::pow(v_plus+v_minus,2)+
             1.0/4.0*std::pow(d_m/A,2)*(del_z[j-1]*d_friction/D + d_K)*(dvdh_plus + dvdh_minus);
          if (d_printJacobian) d_Jacobian[j][0] = J_j_0;

          // compute Jacobian entry J[j][j]
          double J_j_j = -1.0*std::pow(d_m/A,2)*dvdp_minus - 2.0*g*del_z[j-1]*std::cos(d_theta)*
             dvdp_minus/std::pow(v_plus+v_minus,2)+
             1.0/4.0*std::pow(d_m/A,2)*(del_z[j-1]*d_friction/D + d_K)*dvdp_minus - 1.0;
          if (d_printJacobian) d_Jacobian[j][j] = J_j_j;

          // compute Jacobian entry J[j][j+1]
          double J_j_jp1 = std::pow(d_m/A,2)*dvdp_plus - 2.0*g*del_z[j-1]*std::cos(d_theta)*
             dvdp_plus/std::pow(v_plus+v_minus,2)+
             1.0/4.0*std::pow(d_m/A,2)*(del_z[j-1]*d_friction/D + d_K)*dvdp_plus + 1.0;
          if (d_printJacobian) d_Jacobian[j][j+1] = J_j_jp1;
        
          // get delta solution variables
          double du_j   = inputVec->getValueByLocalID(j);   // du[j]
          double du_jp1 = inputVec->getValueByLocalID(j+1); // du[j+1]

          // compute entry j of Jacobian matrix-vector product J*du
          double R_j = J_j_0*du_0 + J_j_j*du_j + J_j_jp1*du_jp1;

          // put residual value in residual vector
          outputVec->setValueByLocalID(j, R_j);
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

// derivative of enthalpy with respect to pressure
double SubchannelTwoEqLinearOperator::dhdp(double T, double p){
      // calculate perturbation
      double b = pow(d_machinePrecision,0.5);
      double pert = (1.0 + p)*b; // perturbation

      // calculate perturbed value
      std::map<std::string, boost::shared_ptr<std::vector<double> > > enthalpyArgMap_pert;
      enthalpyArgMap_pert.insert(std::make_pair("temperature",new std::vector<double>(1,T)));
      enthalpyArgMap_pert.insert(std::make_pair("pressure",   new std::vector<double>(1,p+pert)));
      std::vector<double> enthalpyResult_pert(1);
      d_subchannelPhysicsModel->getProperty("Enthalpy",enthalpyResult_pert,enthalpyArgMap_pert); 
      double h_pert = enthalpyResult_pert[0];

      // calculate unperturbed value
      std::map<std::string, boost::shared_ptr<std::vector<double> > > enthalpyArgMap;
      enthalpyArgMap.insert(std::make_pair("temperature",new std::vector<double>(1,T)));
      enthalpyArgMap.insert(std::make_pair("pressure",   new std::vector<double>(1,p)));
      std::vector<double> enthalpyResult(1);
      d_subchannelPhysicsModel->getProperty("Enthalpy",enthalpyResult,enthalpyArgMap); 
      double h = enthalpyResult[0];

      // calculate derivative
      return (h_pert - h)/pert;
}

// derivative of specific volume with respect to enthalpy
double SubchannelTwoEqLinearOperator::dvdh(double h, double p){
      // calculate perturbation
      double b = pow(d_machinePrecision,0.5);
      double pert = (1.0 + h)*b; // perturbation

      // calculate perturbed value
      std::map<std::string, boost::shared_ptr<std::vector<double> > > specificVolumeArgMap_pert;
      specificVolumeArgMap_pert.insert(std::make_pair("enthalpy",new std::vector<double>(1,h+pert)));
      specificVolumeArgMap_pert.insert(std::make_pair("pressure",new std::vector<double>(1,p)));
      std::vector<double> specificVolumeResult_pert(1);
      d_subchannelPhysicsModel->getProperty("SpecificVolume",specificVolumeResult_pert,specificVolumeArgMap_pert); 
      double v_pert = specificVolumeResult_pert[0];

      // calculate unperturbed value
      std::map<std::string, boost::shared_ptr<std::vector<double> > > specificVolumeArgMap;
      specificVolumeArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,h)));
      specificVolumeArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,p)));
      std::vector<double> specificVolumeResult(1);
      d_subchannelPhysicsModel->getProperty("SpecificVolume",specificVolumeResult,specificVolumeArgMap); 
      double v = specificVolumeResult[0];

      // calculate derivative
      return (v_pert - v)/pert;
}

// derivative of specific volume with respect to pressure
double SubchannelTwoEqLinearOperator::dvdp(double h, double p){
      // calculate perturbation
      double b = pow(d_machinePrecision,0.5);
      double pert = (1.0 + p)*b; // perturbation

      // calculate perturbed value
      std::map<std::string, boost::shared_ptr<std::vector<double> > > specificVolumeArgMap_pert;
      specificVolumeArgMap_pert.insert(std::make_pair("enthalpy",new std::vector<double>(1,h)));
      specificVolumeArgMap_pert.insert(std::make_pair("pressure",new std::vector<double>(1,p+pert)));
      std::vector<double> specificVolumeResult_pert(1);
      d_subchannelPhysicsModel->getProperty("SpecificVolume",specificVolumeResult_pert,specificVolumeArgMap_pert); 
      double v_pert = specificVolumeResult_pert[0];

      // calculate unperturbed value
      std::map<std::string, boost::shared_ptr<std::vector<double> > > specificVolumeArgMap;
      specificVolumeArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,h)));
      specificVolumeArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,p)));
      std::vector<double> specificVolumeResult(1);
      d_subchannelPhysicsModel->getProperty("SpecificVolume",specificVolumeResult,specificVolumeArgMap); 
      double v = specificVolumeResult[0];

      // calculate derivative
      return (v_pert - v)/pert;
}

AMP::LinearAlgebra::Vector::shared_ptr  SubchannelTwoEqLinearOperator::subsetInputVector(AMP::LinearAlgebra::Vector::shared_ptr vec)
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

AMP::LinearAlgebra::Vector::shared_ptr  SubchannelTwoEqLinearOperator::subsetOutputVector(AMP::LinearAlgebra::Vector::shared_ptr vec)
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

}
}
