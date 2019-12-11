
#include "ElementPhysicsModelFactory.h"

#ifdef USE_EXT_LIBMESH
// mechanics material models
#include "AMP/operators/mechanics/ElasticDamageThermalStrainModel.h"
#include "AMP/operators/mechanics/GeneralCladThermalCreepPlasticModel.h"
#include "AMP/operators/mechanics/IsotropicElasticModel.h"
#include "AMP/operators/mechanics/PericElastoViscoPlasticModel.h"
#include "AMP/operators/mechanics/ThermalStrainMaterialModel.h"
#include "AMP/operators/mechanics/ThermalVonMisesMatModel.h"
#include "AMP/operators/mechanics/VonMisesElastoPlasticModel.h"
#include "AMP/operators/mechanics/VonMises_IsotropicKinematicHardening.h"

// flow transport model
#include "AMP/operators/flow/FlowTransportModel.h"

// diffusion transport model
#include "AMP/operators/diffusion/DiffusionTransportModel.h"

// diffusion transport tensor model
#include "AMP/operators/diffusion/DiffusionTransportTensorModel.h"

// diffusion transport cylindrical model
#include "AMP/operators/diffusion/DiffusionCylindricalTransportModel.h"

// Pellet Contact Conductance model
#include "AMP/operators/libmesh/PelletContactConductanceModel.h"

// Convective Heat Coefficient model
#include "AMP/operators/subchannel/ConvectiveHeatCoefficient.h"

// source physics model
#include "AMP/operators/libmesh/SourcePhysicsModel.h"

// mass density model
#include "AMP/operators/libmesh/MassDensityModel.h"

// manufactured diffusion transport model
#include "AMP/operators/ManufacturedDiffusionTransportModel.h"

// subchannel physics model
#include "AMP/operators/subchannel/SubchannelPhysicsModel.h"
#endif


#define resetElementPhysicsModel( NAME )                        \
    do {                                                        \
        if ( name == #NAME )                                    \
            retElementPhysicsModel.reset( new NAME( params ) ); \
    } while ( 0 )


namespace AMP {
namespace Operator {


std::shared_ptr<ElementPhysicsModel> ElementPhysicsModelFactory::createElementPhysicsModel(
    std::shared_ptr<Database> elementPhysicsModelDb )
{
    std::shared_ptr<ElementPhysicsModel> retElementPhysicsModel;
    std::shared_ptr<ElementPhysicsModelParameters> params;

    AMP_INSIST(
        elementPhysicsModelDb.get() != nullptr,
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

    AMP_INSIST( retElementPhysicsModel != nullptr, "requested model " + name + " is invalid" );
    return retElementPhysicsModel;
}


} // namespace Operator
} // namespace AMP
