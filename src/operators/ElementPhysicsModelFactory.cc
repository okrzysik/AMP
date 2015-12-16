
#include "ElementPhysicsModelFactory.h"

#ifdef USE_EXT_LIBMESH
// mechanics material models
#include "operators/mechanics/ElasticDamageThermalStrainModel.h"
#include "operators/mechanics/GeneralCladThermalCreepPlasticModel.h"
#include "operators/mechanics/IsotropicElasticModel.h"
#include "operators/mechanics/PericElastoViscoPlasticModel.h"
#include "operators/mechanics/ThermalStrainMaterialModel.h"
#include "operators/mechanics/ThermalVonMisesMatModel.h"
#include "operators/mechanics/VonMisesElastoPlasticModel.h"
#include "operators/mechanics/VonMises_IsotropicKinematicHardening.h"

// flow transport model
#include "operators/flow/FlowTransportModel.h"

// diffusion transport model
#include "operators/diffusion/DiffusionTransportModel.h"

// diffusion transport tensor model
#include "operators/diffusion/DiffusionTransportTensorModel.h"

// diffusion transport cylindrical model
#include "operators/diffusion/DiffusionCylindricalTransportModel.h"

// Pellet Contact Conductance model
#include "operators/libmesh/PelletContactConductanceModel.h"

// Convective Heat Coefficient model
#include "operators/subchannel/ConvectiveHeatCoefficient.h"

// source physics model
#include "operators/libmesh/SourcePhysicsModel.h"

// mass density model
#include "operators/libmesh/MassDensityModel.h"

// manufactured diffusion transport model
#include "operators/ManufacturedDiffusionTransportModel.h"

// subchannel physics model
#include "operators/subchannel/SubchannelPhysicsModel.h"
#endif


#define resetElementPhysicsModel( NAME )                                         \
    do {                                                                         \
        if ( name == #NAME ) retElementPhysicsModel.reset( new NAME( params ) ); \
    } while ( 0 )


namespace AMP {
namespace Operator {


AMP::shared_ptr<ElementPhysicsModel> ElementPhysicsModelFactory::createElementPhysicsModel(
    AMP::shared_ptr<Database> elementPhysicsModelDb )
{
    AMP::shared_ptr<ElementPhysicsModel> retElementPhysicsModel;
    AMP::shared_ptr<ElementPhysicsModelParameters> params;

    AMP_INSIST(
        elementPhysicsModelDb.get() != NULL,
        "ElementPhysicsModelFactory::createElementPhysicsModel:: NULL Database object input" );

    std::string name = elementPhysicsModelDb->getString( "name" );

    params.reset( new ElementPhysicsModelParameters( elementPhysicsModelDb ) );

#ifdef USE_EXT_LIBMESH
    resetElementPhysicsModel( IsotropicElasticModel );
    resetElementPhysicsModel( ThermalStrainMaterialModel );
    resetElementPhysicsModel( VonMisesElastoPlasticModel );
    resetElementPhysicsModel( ThermalVonMisesMatModel );
    resetElementPhysicsModel( GeneralCladThermalCreepPlasticModel );
    resetElementPhysicsModel( VonMises_IsotropicKinematicHardening );
    resetElementPhysicsModel( ElasticDamageThermalStrainModel );
    resetElementPhysicsModel( PericElastoViscoPlasticModel );
    resetElementPhysicsModel( DiffusionTransportModel );
    resetElementPhysicsModel( DiffusionTransportTensorModel );
    resetElementPhysicsModel( DiffusionCylindricalTransportModel );
    resetElementPhysicsModel( FlowTransportModel );
    resetElementPhysicsModel( PelletContactConductanceModel );
    resetElementPhysicsModel( ConvectiveHeatCoefficient );
    resetElementPhysicsModel( SourcePhysicsModel );
    resetElementPhysicsModel( MassDensityModel );
    resetElementPhysicsModel( ManufacturedDiffusionTransportModel );
    resetElementPhysicsModel( SubchannelPhysicsModel );
#endif

    AMP_INSIST( retElementPhysicsModel != NULL, "requested model " + name + " is invalid" );
    return retElementPhysicsModel;
}


} // namespace Operator
} // namespace AMP
