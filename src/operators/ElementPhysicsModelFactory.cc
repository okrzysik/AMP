#include "ElementPhysicsModelFactory.h"

// mechanics material models
#include "IsotropicElasticModel.h"
#include "ThermalStrainMaterialModel.h"
#include "VonMisesElastoPlasticModel.h"
#include "ThermalVonMisesMatModel.h"
#include "GeneralCladThermalCreepPlasticModel.h"
#include "VonMises_IsotropicKinematicHardening.h"
#include "ElasticDamageThermalStrainModel.h"
#include "PericElastoViscoPlasticModel.h"

// flow transport model
#include "FlowTransportModel.h"

// diffusion transport model
#include "DiffusionTransportModel.h"

// diffusion transport tensor model
#include "DiffusionTransportTensorModel.h"

// Pellet Contact Conductance model
#include "PelletContactConductanceModel.h"

// source physics model
#include "SourcePhysicsModel.h"

// mass density model
#include "MassDensityModel.h"

// manufactured diffusion transport model
#include "ManufacturedDiffusionTransportModel.h"

namespace AMP {
namespace Operator {

boost::shared_ptr<ElementPhysicsModel>
ElementPhysicsModelFactory::createElementPhysicsModel(boost::shared_ptr<Database>  elementPhysicsModelDb)
{
    boost::shared_ptr<ElementPhysicsModel> retElementPhysicsModel;
    boost::shared_ptr<ElementPhysicsModelParameters> params;
    
    AMP_INSIST(elementPhysicsModelDb.get()!=NULL, "ElementPhysicsModelFactory::createElementPhysicsModel:: NULL Database object input");

    std::string name = elementPhysicsModelDb->getString("name");
    
    params.reset( new ElementPhysicsModelParameters(elementPhysicsModelDb));

    if(name=="IsotropicElasticModel")
      {
    retElementPhysicsModel.reset(new IsotropicElasticModel(params));
      }
    else if(name=="ThermalStrainMaterialModel")
      {
    retElementPhysicsModel.reset(new ThermalStrainMaterialModel(params));
      }
    else if(name=="VonMisesElastoPlasticModel")
      {
    retElementPhysicsModel.reset(new VonMisesElastoPlasticModel(params));
      }
    else if(name=="ThermalVonMisesMatModel")
      {
    retElementPhysicsModel.reset(new ThermalVonMisesMatModel(params));
      }
    else if(name=="GeneralCladThermalCreepPlasticModel")
      {
        retElementPhysicsModel.reset(new GeneralCladThermalCreepPlasticModel(params));
      }
    else if(name=="VonMises_IsotropicKinematicHardening")
      {
        retElementPhysicsModel.reset(new VonMises_IsotropicKinematicHardening(params));
      }
    else if(name=="ElasticDamageThermalStrainModel")
      {
        retElementPhysicsModel.reset(new ElasticDamageThermalStrainModel(params));
      }
    else if(name=="PericElastoViscoPlasticModel")
      {
        retElementPhysicsModel.reset(new PericElastoViscoPlasticModel(params));
      }
    else if(name=="DiffusionTransportModel")
      {
        retElementPhysicsModel.reset(new DiffusionTransportModel(params));
      }
    else if(name=="DiffusionTransportTensorModel")
      {
        retElementPhysicsModel.reset(new DiffusionTransportTensorModel(params));
      }
    else if(name=="FlowTransportModel")
      {
        retElementPhysicsModel.reset(new FlowTransportModel(params));
      }
    else if(name=="PelletContactConductanceModel")
      {
    retElementPhysicsModel.reset(new PelletContactConductanceModel(params));
      }
    else if(name=="SourcePhysicsModel")
      {
    retElementPhysicsModel.reset(new SourcePhysicsModel(params));
      }
    else if(name=="MassDensityModel")
      {
    retElementPhysicsModel.reset(new MassDensityModel(params));
      }
    else if(name=="ManufacturedDiffusionTransportModel")
      {
    retElementPhysicsModel.reset(new ManufacturedDiffusionTransportModel(params));
      } 
    else
    {
        AMP_INSIST(false, "requested model "+name+" is invalid");
    }

    return retElementPhysicsModel;

}
  
}
}

