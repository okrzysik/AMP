
#include "ElementPhysicsModelFactory.h"

// mechanics material models
#include "mechanics/IsotropicElasticModel.h"
#include "mechanics/ThermalStrainMaterialModel.h"
#include "mechanics/VonMisesElastoPlasticModel.h"
#include "mechanics/ThermalVonMisesMatModel.h"
#include "mechanics/GeneralCladThermalCreepPlasticModel.h"
#include "mechanics/VonMises_IsotropicKinematicHardening.h"
#include "mechanics/ElasticDamageThermalStrainModel.h"
#include "mechanics/PericElastoViscoPlasticModel.h"

// flow transport model
#include "flow/FlowTransportModel.h"

// diffusion transport model
#include "diffusion/DiffusionTransportModel.h"

// diffusion transport tensor model
#include "diffusion/DiffusionTransportTensorModel.h"

// diffusion transport cylindrical model
#include "diffusion/DiffusionCylindricalTransportModel.h"

// Pellet Contact Conductance model
#include "PelletContactConductanceModel.h"

// source physics model
#include "SourcePhysicsModel.h"

// mass density model
#include "MassDensityModel.h"

// manufactured diffusion transport model
#include "ManufacturedDiffusionTransportModel.h"

// subchannel physics model
#include "SubchannelPhysicsModel.h"

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
        else if(name=="DiffusionCylindricalTransportModel")
        {
          retElementPhysicsModel.reset(new DiffusionCylindricalTransportModel(params));
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
        else if(name=="SubchannelPhysicsModel")
        {
          retElementPhysicsModel.reset(new SubchannelPhysicsModel(params));
        } 
        else
        {
          AMP_INSIST(false, "requested model "+name+" is invalid");
        }

        return retElementPhysicsModel;

      }

  }
}

