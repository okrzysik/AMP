#include <cstring>

#include "SourcePhysicsModel.h"
#include "ElementPhysicsModel.h"
#include "materials/Material.h"
#include "ElementPhysicsModelParameters.h"
//#include "diffusion/DiffusionTransportModel.h"
#include "MassDensityModel.h"
#include "boost/shared_ptr.hpp"

#include "ManufacturedSourceModel1.h"
#include "ManufacturedSourceModel2.h"


/* Libmesh files */
#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
#include "quadrature_gauss.h"


namespace AMP {
namespace Operator {


SourcePhysicsModel::SourcePhysicsModel (const boost::shared_ptr<SourcePhysicsModelParameters>& params )
        : ElementPhysicsModel(params){
         d_useMaterialsLibrary = (params->d_db)->getBoolWithDefault("USE_MATERIALS_LIBRARY",false);

         if(d_useMaterialsLibrary == true) 
         {
             AMP_INSIST( (params->d_db->keyExists("Material")), "Key ''Material'' is missing!" );
             std::string matname = params->d_db->getString("Material");
             d_material = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create(matname);
         }

         elementPhysicsParams.reset(new AMP::Operator::ElementPhysicsModelParameters( params->d_db ));
         d_physicsName = elementPhysicsParams->d_db->getString("USE_ELEMENT_PHYSICS");

         d_DefaultTemperature = params->d_db->getDoubleWithDefault("Default_Temperature",300.);
         d_DefaultConcentration = params->d_db->getDoubleWithDefault("Default_Concentration",0.1);
         d_DefaultBurnup = params->d_db->getDoubleWithDefault("Default_Burnup",0.);
         
         d_defaults.resize(3);

         d_defaults[0] = d_DefaultTemperature;
         d_defaults[1] = d_DefaultConcentration;
         d_defaults[2] = d_DefaultBurnup;

    }


void SourcePhysicsModel::getConstitutiveProperty(       std::vector<double>                & result, 
                                         std::vector< std::vector<double> > & InputVec, 
                                         std::vector< std::vector<double> > & ,
                                      const std::vector<Point> & Coordinates                )
     {
AMP_ERROR("SourcePhysicsModel is not converted yet");
/*
         if(d_physicsName=="DiffusionTransportModel")
         {
             elementPhysicsModel.reset(new AMP::Operator::DiffusionTransportModel(boost::dynamic_pointer_cast<DiffusionTransportModelParameters>(elementPhysicsParams) )); 
             boost::shared_ptr<DiffusionTransportModel> tmp = boost::dynamic_pointer_cast<DiffusionTransportModel>(elementPhysicsModel);  
             d_property = tmp->getProperty(); 

             std::map<std::string, boost::shared_ptr<std::vector<double> > > inputMaterialParameters;

             std::string temperatureString = "temperature"; // in the future get from input file
             std::string burnupString = "burnup"; // in the future get from input file
             std::string oxygenString = "concentration"; // in the future get from input file

             boost::shared_ptr<std::vector<double> > tempVec(new std::vector<double> );      
             boost::shared_ptr<std::vector<double> > burnupVec( new std::vector<double>(result.size(),d_DefaultBurnup) );
             boost::shared_ptr<std::vector<double> > oxygenVec( new std::vector<double>(result.size(),d_DefaultConcentration) );

             inputMaterialParameters.insert( std::make_pair( temperatureString, tempVec) );
             inputMaterialParameters.insert( std::make_pair( burnupString, burnupVec) );
             inputMaterialParameters.insert( std::make_pair( oxygenString, oxygenVec) );

             int switchId = InputVec.size();  
             switch(switchId) {
                 case 1 : {
                          (*tempVec) = InputVec[0];
                          tmp->getTransport(result, inputMaterialParameters, Coordinates);
                          return;
                      }
                 default : {
                          std::vector<double> temp(result.size(),d_DefaultTemperature);
                          (*tempVec) = temp;
                          tmp->getTransport(result, inputMaterialParameters, Coordinates);
                          return;
                      }
             }

         }
         else if (d_physicsName=="MassDensityModel")
         {
             elementPhysicsModel.reset(new AMP::Operator::MassDensityModel(boost::dynamic_pointer_cast<MassDensityModelParameters>(elementPhysicsParams) )); 
             boost::shared_ptr<MassDensityModel> tmp = boost::dynamic_pointer_cast<MassDensityModel>(elementPhysicsModel);  

             std::string eqnname = elementPhysicsParams->d_db->getString("Equation");
             AMP_INSIST((eqnname == "ThermalSource"), "SourcePhysicsModel should be implemented by User for this Equation"); 

             int switchId = InputVec.size();  

             switch(switchId) {
                 case 1 : {
                          std::vector<double> DefaultVec0(result.size());
                          std::vector<double> DefaultVec1(result.size());
                          std::vector<double> DefaultVec2(result.size());
                          for (size_t i=0; i<result.size(); i++){ 
                              DefaultVec0[i]=d_defaults[0];
                              DefaultVec1[i]=d_defaults[1];
                              DefaultVec2[i]=d_defaults[2];
                          }
                          tmp->getDensityMechanics(result, DefaultVec0, DefaultVec1, DefaultVec2);
                          for (size_t i=0; i<result.size(); i++){
                              result[i] = result[i]*InputVec[0][i];
                          } 
                          break;
                      }
                 case 2 : {
                          std::vector<double> DefaultVec(result.size());
                          for (size_t i=0; i<result.size(); i++) {DefaultVec[i]=d_defaults[2];}
                          tmp->getDensityMechanics(result, InputVec[0], InputVec[1], DefaultVec);
                          break;
                      }
                 case 3 : {
                          tmp->getDensityMechanics(result, InputVec[0], InputVec[1], InputVec[2]);
                          break;
                      }
                 default:
                      assert(false);
             }
         }
         else if (d_physicsName=="ManufacturedSourceModel1")
            {
                elementPhysicsModel.reset(new AMP::Operator::ManufacturedSourceModel1(boost::dynamic_pointer_cast<ManufacturedSourceModel1Parameters>(elementPhysicsParams) )); 
                boost::shared_ptr<ManufacturedSourceModel1> tmp = boost::dynamic_pointer_cast<ManufacturedSourceModel1>(elementPhysicsModel);  
                
                int switchId = InputVec.size();  
                switch(switchId) {
                    case 1 : {
                        std::vector<double> DefaultVec;

                        DefaultVec.resize(result.size());
                        for (size_t i=0; i<result.size(); i++){ 
                            DefaultVec[i] = InputVec[0][i];
                        }
                        tmp->getManufacturedSource1(result, DefaultVec, Coordinates);
                        break;
                    }
                    default:
                        assert(false);
                }
                
            }
            else if (d_physicsName=="ManufacturedSourceModel2")
            {
                elementPhysicsModel.reset(new AMP::Operator::ManufacturedSourceModel2(boost::dynamic_pointer_cast<ManufacturedSourceModel2Parameters>(elementPhysicsParams) )); 
                boost::shared_ptr<ManufacturedSourceModel2> tmp = boost::dynamic_pointer_cast<ManufacturedSourceModel2>(elementPhysicsModel);  
                
                int switchId = InputVec.size();  
                switch(switchId) {
                    case 1 : {
                        std::vector<double> DefaultVec;
                        
                        DefaultVec.resize(result.size());
                        for (size_t i=0; i<result.size(); i++){ 
                            DefaultVec[i] = InputVec[0][i];
                        }
                        tmp->getManufacturedSource2(result, DefaultVec, Coordinates);
                        break;
                    }
                    default:
                        assert(false);
                }
                
            }
*/
}


}
}


