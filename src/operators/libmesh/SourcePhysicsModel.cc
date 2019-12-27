#include "AMP/operators/libmesh/SourcePhysicsModel.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/ElementPhysicsModel.h"
#include "AMP/operators/ElementPhysicsModelParameters.h"
#include "AMP/operators/ManufacturedSourceModel1.h"
#include "AMP/operators/ManufacturedSourceModel2.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/MassDensityModel.h"
#include <memory>

#include "ProfilerApp.h"

#include <cstring>


namespace AMP {
namespace Operator {


SourcePhysicsModel::SourcePhysicsModel(
    const std::shared_ptr<SourcePhysicsModelParameters> &params )
    : ElementPhysicsModel( params )
{
    d_useMaterialsLibrary = ( params->d_db )->getWithDefault( "USE_MATERIALS_LIBRARY", false );

    if ( d_useMaterialsLibrary == true ) {
        AMP_INSIST( ( params->d_db->keyExists( "Material" ) ), "Key ''Material'' is missing!" );
        std::string matname = params->d_db->getString( "Material" );
        d_material = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( matname );
    } else {
        // Read constant value for the associated property.
        d_constantProperty = params->d_db->getWithDefault<double>( "CONSTANT_VALUE", 1.0 );
    }

    d_elementPhysicsParams.reset(
        new AMP::Operator::ElementPhysicsModelParameters( params->d_db ) );
    d_physicsName = d_elementPhysicsParams->d_db->getString( "USE_ELEMENT_PHYSICS" );

    d_DefaultTemperature   = params->d_db->getWithDefault<double>( "Default_Temperature", 300. );
    d_DefaultConcentration = params->d_db->getWithDefault<double>( "Default_Concentration", 0.1 );
    d_DefaultBurnup        = params->d_db->getWithDefault<double>( "Default_Burnup", 0. );

    d_defaults.resize( 3 );

    d_defaults[0] = d_DefaultTemperature;
    d_defaults[1] = d_DefaultConcentration;
    d_defaults[2] = d_DefaultBurnup;

    if ( d_physicsName == "DiffusionTransportModel" ) {
        d_elementPhysicsModel.reset( new AMP::Operator::DiffusionTransportModel(
            std::dynamic_pointer_cast<DiffusionTransportModelParameters>(
                d_elementPhysicsParams ) ) );
        std::shared_ptr<DiffusionTransportModel> tmp =
            std::dynamic_pointer_cast<DiffusionTransportModel>( d_elementPhysicsModel );
        d_property = tmp->getProperty();

        const std::string temperatureString = "temperature";   // in the future get from input file
        const std::string burnupString      = "burnup";        // in the future get from input file
        const std::string oxygenString      = "concentration"; // in the future get from input file

        std::shared_ptr<std::vector<double>> tempVec(
            new std::vector<double>( 1, d_DefaultTemperature ) );
        std::shared_ptr<std::vector<double>> burnupVec(
            new std::vector<double>( 1, d_DefaultBurnup ) );
        std::shared_ptr<std::vector<double>> oxygenVec(
            new std::vector<double>( 1, d_DefaultConcentration ) );

        d_inputMaterialParameters.clear();
        d_inputMaterialParameters.insert( std::make_pair( temperatureString, tempVec ) );
        d_inputMaterialParameters.insert( std::make_pair( burnupString, burnupVec ) );
        d_inputMaterialParameters.insert( std::make_pair( oxygenString, oxygenVec ) );
    } else if ( d_physicsName == "MassDensityModel" ) {
        if ( d_useMaterialsLibrary ) {
            d_elementPhysicsModel.reset( new AMP::Operator::MassDensityModel(
                std::dynamic_pointer_cast<MassDensityModelParameters>( d_elementPhysicsParams ) ) );
        }
    } else if ( d_physicsName == "ManufacturedSourceModel1" ) {
        d_elementPhysicsModel.reset( new AMP::Operator::ManufacturedSourceModel1(
            std::dynamic_pointer_cast<ManufacturedSourceModel1Parameters>(
                d_elementPhysicsParams ) ) );
    } else if ( d_physicsName == "ManufacturedSourceModel2" ) {
        d_elementPhysicsModel.reset( new AMP::Operator::ManufacturedSourceModel2(
            std::dynamic_pointer_cast<ManufacturedSourceModel2Parameters>(
                d_elementPhysicsParams ) ) );
    }
}


void SourcePhysicsModel::getConstitutiveProperty( std::vector<double> &result,
                                                  const std::vector<std::vector<double>> &InputVec,
                                                  const std::vector<std::vector<double>> &,
                                                  const std::vector<libMesh::Point> &Coordinates )
{
    PROFILE_START( "getConstitutiveProperty", 6 );
    if ( d_physicsName == "DiffusionTransportModel" ) {
        std::shared_ptr<DiffusionTransportModel> tmp =
            std::dynamic_pointer_cast<DiffusionTransportModel>( d_elementPhysicsModel );

        auto it = d_inputMaterialParameters.find( "temperature" );
        AMP_ASSERT( it != d_inputMaterialParameters.end() );
        std::vector<double> &tempVec = *( it->second );

        if ( tempVec.size() != result.size() ) {
            // Resize the vectors to match the number of outputs
            std::vector<double> &burnupVec =
                *( d_inputMaterialParameters.find( "burnup" )->second );
            std::vector<double> &oxygenVec =
                *( d_inputMaterialParameters.find( "concentration" )->second );
            tempVec.resize( result.size(), d_DefaultTemperature );
            burnupVec.resize( result.size(), d_DefaultBurnup );
            oxygenVec.resize( result.size(), d_DefaultConcentration );
        }
        if ( InputVec.size() == 1 ) {
            AMP_ASSERT( InputVec[0].size() == result.size() );
            for ( size_t i = 0; i < tempVec.size(); i++ )
                tempVec[i] = InputVec[0][i];
        } else {
            for ( auto &elem : tempVec )
                elem = d_DefaultTemperature;
        }

        tmp->getTransport( result, d_inputMaterialParameters, Coordinates );
    } else if ( d_physicsName == "MassDensityModel" ) {
        if ( !d_useMaterialsLibrary ) {
            for ( size_t i = 0; i < result.size(); i++ ) {
                result[i] = d_constantProperty * InputVec[0][i];
            }
        } else {
            std::shared_ptr<MassDensityModel> tmp =
                std::dynamic_pointer_cast<MassDensityModel>( d_elementPhysicsModel );

            std::string eqnname = d_elementPhysicsParams->d_db->getString( "Equation" );
            AMP_INSIST( ( eqnname == "ThermalSource" ),
                        "SourcePhysicsModel should be implemented by User for this Equation" );

            if ( InputVec.size() == 1 ) {
                std::vector<double> DefaultVec0( result.size() );
                std::vector<double> DefaultVec1( result.size() );
                std::vector<double> DefaultVec2( result.size() );
                for ( size_t i = 0; i < result.size(); i++ ) {
                    DefaultVec0[i] = d_defaults[0];
                    DefaultVec1[i] = d_defaults[1];
                    DefaultVec2[i] = d_defaults[2];
                }
                tmp->getDensityMechanics( result, DefaultVec0, DefaultVec1, DefaultVec2 );
                for ( size_t i = 0; i < result.size(); i++ ) {
                    result[i] = result[i] * InputVec[0][i];
                }
            } else {
                AMP_ERROR( "Invalid case" );
            }
        }
        // AMP::pout<<" Power [W/m^3]: "<< result[0] << " Specific Power[W/kg]: "<<InputVec[0][0]
        // <<" Density [kg/m^3]:
        // " << result[0]/InputVec[0][0] <<std::endl;
    } else if ( d_physicsName == "ManufacturedSourceModel1" ) {
        std::shared_ptr<ManufacturedSourceModel1> tmp =
            std::dynamic_pointer_cast<ManufacturedSourceModel1>( d_elementPhysicsModel );

        if ( InputVec.size() == 1 ) {
            std::vector<double> DefaultVec;

            DefaultVec.resize( result.size() );
            for ( size_t i = 0; i < result.size(); i++ ) {
                DefaultVec[i] = InputVec[0][i];
            }
            tmp->getManufacturedSource1( result, DefaultVec, Coordinates );
        } else {
            AMP_ERROR( "Invalid case" );
        }
    } else if ( d_physicsName == "ManufacturedSourceModel2" ) {
        std::shared_ptr<ManufacturedSourceModel2> tmp =
            std::dynamic_pointer_cast<ManufacturedSourceModel2>( d_elementPhysicsModel );

        if ( InputVec.size() == 1 ) {
            std::vector<double> DefaultVec;

            DefaultVec.resize( result.size() );
            for ( size_t i = 0; i < result.size(); i++ ) {
                DefaultVec[i] = InputVec[0][i];
            }
            tmp->getManufacturedSource2( result, DefaultVec, Coordinates );
        } else {
            AMP_ERROR( "Invalid case" );
        }
    }
    PROFILE_STOP( "getConstitutiveProperty", 6 );
}
} // namespace Operator
} // namespace AMP
