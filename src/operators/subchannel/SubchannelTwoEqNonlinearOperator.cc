
#include "operators/subchannel/SubchannelTwoEqNonlinearOperator.h"
#include "operators/subchannel/SubchannelOperatorParameters.h"
#include "ampmesh/StructuredMeshHelper.h"
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

      d_params = myparams;

      d_Pout = getDoubleParameter(myparams,"Exit_Pressure",15.5132e6);
      d_Tin  = getDoubleParameter(myparams,"Inlet_Temperature",569.26);  
      d_m    = getDoubleParameter(myparams,"Mass_Flow_Rate",0.3522);  
      d_gamma     = getDoubleParameter(myparams,"Fission_Heating_Coefficient",0.0);  
      d_theta     = getDoubleParameter(myparams,"Channel_Angle",0.0);  
      d_channelDia= getDoubleParameter(myparams,"Channel_Diameter",0.0);  
      d_reynolds  = getDoubleParameter(myparams,"Reynolds",0.0);  
      d_prandtl   = getDoubleParameter(myparams,"Prandtl",0.0);  
      d_friction  = getDoubleParameter(myparams,"Friction_Factor",0.1);  
      d_pitch     = getDoubleParameter(myparams,"Lattice_Pitch",0.0128016);  
      d_diameter  = getDoubleParameter(myparams,"Rod_Diameter",0.0097028);  
      d_K    = getDoubleParameter(myparams,"Form_Loss_Coefficient",0.2);  
      d_source = getStringParameter(myparams,"Heat_Source_Type","totalHeatGeneration");

      // get additional parameters based on heat source type
      if (d_source == "totalHeatGeneration") {
          d_Q    = getDoubleParameter(myparams,"Rod_Power",66.81e3);  
          d_heatShape = getStringParameter(myparams,"Heat_Shape","Sinusoidal");
      }else if (d_source == "averageCladdingTemperature"){
        d_channelFractions = (myparams->d_db)->getDoubleArray("ChannelFractions");
      }

      d_subchannelPhysicsModel   = myparams->d_subchannelPhysicsModel;  
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

void SubchannelTwoEqNonlinearOperator::fillSubchannelGrid(AMP::Mesh::Mesh::shared_ptr mesh)
{
    // Create the grid for all processors
    std::set<double> x, y, z;
    if ( mesh.get() != NULL ) {
        AMP::Mesh::MeshIterator it = mesh->getIterator( AMP::Mesh::Vertex, 0 );
        for (size_t i=0; i<it.size(); i++) {
            std::vector<double> coord = it->coord();
            AMP_ASSERT(coord.size()==3);
            x.insert( coord[0] );
            y.insert( coord[1] );
            z.insert( coord[2] );
            ++it;
        }
    }
    d_Mesh->getComm().setGather(x);
    d_Mesh->getComm().setGather(y);
    d_Mesh->getComm().setGather(z);
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
    d_x = std::vector<double>(x.begin(),x.end());
    d_y = std::vector<double>(y.begin(),y.end());
    d_z = std::vector<double>(z.begin(),z.end());
    size_t Nx = d_x.size()-1;
    size_t Ny = d_y.size()-1;
    size_t Nz = d_z.size()-1;
    if ( mesh.get() != NULL ) 
        AMP_ASSERT(Nx*Ny*Nz==mesh->numGlobalElements(AMP::Mesh::Volume));
    d_numSubchannels = Nx*Ny;
}

int SubchannelTwoEqNonlinearOperator::getSubchannelIndex( double x, double y )
{
    size_t i = Utilities::findfirst(d_x,x);
    size_t j = Utilities::findfirst(d_y,y);
    if ( i>0 && i<d_x.size() && j>0 && j<d_y.size() )
        return (i-1)+(j-1)*(d_x.size()-1);
    return -1;
}

// apply
void SubchannelTwoEqNonlinearOperator :: apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
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
      const double A = std::pow(d_pitch,2) - pi*std::pow(d_diameter,2)/4.0; // flow area
      const double D = 4.0*A/perimeter;                                     // hydraulic diameter

      // Subset the vectors
      AMP::LinearAlgebra::Vector::const_shared_ptr inputVec = subsetInputVector( u );
      AMP::LinearAlgebra::Vector::shared_ptr outputVec = subsetOutputVector( r );

      AMP::Discretization::DOFManager::shared_ptr dof_manager = inputVec->getDOFManager();
      AMP::Discretization::DOFManager::shared_ptr cladDofManager;
      if (d_source == "averageCladdingTemperature"){
        cladDofManager = d_cladTemperature->getDOFManager();
      }

      fillSubchannelGrid(d_Mesh);

      AMP::Mesh::MeshIterator el = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
      AMP::Mesh::MeshIterator end_el = el.end();

      std::vector<std::vector<AMP::Mesh::MeshElement> > d_elem(d_numSubchannels);
      std::vector<bool> d_ownSubChannel(d_numSubchannels);

      for( ; el != end_el; ++el) {
        std::vector<double> center = el->centroid();
        int index = getSubchannelIndex( center[0], center[1] );
        if ( index>=0 ){
          d_ownSubChannel[index] = true;
          d_elem[index].push_back( *el );
        }
      }//end for el

      for(int isub =0; isub<d_numSubchannels; ++isub){
        if(d_ownSubChannel[isub]){
          boost::shared_ptr<std::vector<AMP::Mesh::MeshElement> > subchannelElements( new std::vector<AMP::Mesh::MeshElement>() );
          subchannelElements->reserve(d_numSubchannels);
          for(size_t ielem=0; ielem<d_elem[isub].size(); ++ielem){
            subchannelElements->push_back(d_elem[isub][ielem]);
          }
          AMP::Mesh::MeshIterator localSubchannelIt = AMP::Mesh::MultiVectorIterator( subchannelElements );
          AMP::Mesh::Mesh::shared_ptr localSubchannel = d_Mesh->Subset( localSubchannelIt  );

          // get the Iterators for the subchannel mesh
          AMP::Mesh::MeshIterator begin_face = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(localSubchannel , 0);
          AMP::Mesh::MeshIterator end_face   = begin_face.end();

          /**
            get boundary values: u has ordering:
            \f[ \vec{u}=\left[\begin{array}{c}
            h_{0^+}\\
            p_{0^+}\\
            \vdots\\
            h_{J^+}\\
            p_{J^+}\\
            \end{array}\right] \f]
            */
          const int numFaces = begin_face.size() ;
          const int numCells = numFaces - 1;

          std::vector<size_t> dofs, scalarDofs;
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
          AMP::Mesh::MeshIterator face = begin_face;
          if (d_source == "averageCladdingTemperature") {
            for (int j=0; j<numCells; j++){
              //evalute flow average temperature 
              dof_manager->getDOFs( face->globalID(), dofs );
              cladDofManager->getDOFs( face->globalID(), scalarDofs );
              double cladMinus = d_cladTemperature->getValueByGlobalID(scalarDofs[0]);
              std::map<std::string, boost::shared_ptr<std::vector<double> > > temperatureArgMap;
              temperatureArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,inputVec->getValueByGlobalID(dofs[0]))));
              temperatureArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,inputVec->getValueByGlobalID(dofs[1]))));
              std::vector<double> flowMinus(1), specificVolMinus(1);
              d_subchannelPhysicsModel->getProperty("Temperature", flowMinus, temperatureArgMap);
              d_subchannelPhysicsModel->getProperty("SpecificVolume", specificVolMinus, temperatureArgMap);

              face++;
              dof_manager->getDOFs( face->globalID(), dofs );
              cladDofManager->getDOFs( face->globalID(), scalarDofs );
              double cladPlus = d_cladTemperature->getValueByGlobalID(scalarDofs[0]);;
              temperatureArgMap.clear();
              temperatureArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,inputVec->getValueByGlobalID(dofs[0]))));
              temperatureArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,inputVec->getValueByGlobalID(dofs[1]))));
              std::vector<double> flowPlus(1), specificVolPlus(1);
              d_subchannelPhysicsModel->getProperty("Temperature", flowPlus, temperatureArgMap); 
              d_subchannelPhysicsModel->getProperty("SpecificVolume", specificVolPlus, temperatureArgMap);

              double cladAvgTemp = (cladMinus + cladPlus)/2.0;
              std::vector<double> flowTemp(1, (flowMinus[0] + flowPlus[0])/2.0);
              std::vector<double> flowDens(1, ((1./specificVolMinus[0]) + (1./specificVolPlus[0]))/2.0);

              std::map<std::string, boost::shared_ptr<std::vector<double> > > convectiveHeatArgMap;
              convectiveHeatArgMap.insert(std::make_pair("temperature",new std::vector<double>(1,flowTemp[0])));
              convectiveHeatArgMap.insert(std::make_pair("density",new std::vector<double>(1,flowDens[0])));
              convectiveHeatArgMap.insert(std::make_pair("diameter",new std::vector<double>(1,d_channelDia)));
              convectiveHeatArgMap.insert(std::make_pair("reynolds",new std::vector<double>(1,d_reynolds)));
              convectiveHeatArgMap.insert(std::make_pair("prandtl",new std::vector<double>(1,d_prandtl)));
              std::vector<double> heff(1); 
              d_subchannelPhysicsModel->getProperty("ConvectiveHeat", heff, convectiveHeatArgMap); 

              double flux = heff[0]*(cladAvgTemp - flowTemp[0]) ;
              double lin  = flux*pi*d_diameter*d_channelFractions[isub] ;
              double flux_sum = 4.0*pi*d_diameter*1.0/4.0*flux;
              double lin_sum = 4.0*d_gamma*1.0/4.0*lin;
              dh[j] = del_z[j] / d_m * (flux_sum + lin_sum);
            }
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
          face = begin_face;
          for(size_t iface = 0; iface < begin_face.size(); ++iface, ++j){

            // ======================================================
            // energy residual
            // ======================================================
            if (face == begin_face){
              /**
                evaluate first residual entry, corresponding to inlet enthalpy:
                \f[ R_0 = h_{in} - h(T_{in},p_{1-})\f]
                */
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
        }
      }//end of isub
      if(f.get() == NULL) {
        outputVec->scale(a);
      } else {
        AMP::LinearAlgebra::Vector::const_shared_ptr fInternal = subsetInputVector( f );
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

  tmp_db->putString("name","SubchannelTwoEqLinearOperator");

  boost::shared_ptr<SubchannelOperatorParameters> outParams(new SubchannelOperatorParameters(tmp_db));
  outParams->d_db = d_params->d_db; 
  outParams->d_dofMap = d_params->d_dofMap; 
  outParams->d_frozenSolution =  subsetInputVector( u );
  outParams->d_subchannelPhysicsModel =  d_subchannelPhysicsModel;

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
