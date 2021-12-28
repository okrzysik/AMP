/*
 * SubchannelPhysicsModel.cc
 *
 */

#include "AMP/operators/subchannel/SubchannelPhysicsModel.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>

namespace AMP::Operator {

SubchannelPhysicsModel::SubchannelPhysicsModel(
    std::shared_ptr<const ElementPhysicsModelParameters> params )
    : ElementPhysicsModel( params )
{
    // get material key
    AMP_INSIST( ( params->d_db->keyExists( "Material" ) ),
                "Subchannel Key ''Material'' is missing!" );
    auto matname = params->d_db->getString( "Material" );
    d_material   = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( matname );

    // get the formulation key
    AMP_INSIST( ( params->d_db->keyExists( "Formulation" ) ),
                "Subchannel Key ''Formulation'' is missing!" );
    auto formulation = params->d_db->getString( "Formulation" );

    // determine which property functions are needed based on formulation key
    std::vector<std::string> properties;
    if ( formulation == std::string( "OneEqnForwardSubstitution" ) ) {
        AMP_ERROR( "The formulation ''OneEqnForwardSubstitution'' has not yet been implemented" );
    } else if ( formulation == std::string( "TwoEqnPicardIteration" ) ) {
        properties.emplace_back( "Enthalpy" );
    } else if ( formulation == std::string( "TwoEqnJFNKPressure" ) ) {
        properties.emplace_back( "Enthalpy" );
        properties.emplace_back( "SpecificVolume" );
    } else if ( formulation == std::string( "TwoEqnJFNKDensity" ) ) {
        properties.emplace_back( "Enthalpy" );
        properties.emplace_back( "Pressure" );
        AMP_ERROR( "The formulation ''TwoEqnJFNKDensity'' has not yet been implemented" );
    } else if ( formulation == std::string( "FunctionsTesting" ) ) {
        properties.emplace_back( "Temperature" );
        properties.emplace_back( "SaturatedLiquidEnthalpy" );
        properties.emplace_back( "SpecificVolume" );
        properties.emplace_back( "ThermalConductivity" );
        properties.emplace_back( "ConvectiveHeat" );
        properties.emplace_back( "DynamicViscosity" );
        properties.emplace_back( "Enthalpy" );
    } else {
        AMP_ERROR( "Invalid Formulation key" );
    }
    // add properties to property pointer map
    for ( auto &propertie : properties ) {
        d_properties.insert( std::make_pair( propertie, d_material->property( propertie ) ) );
    }

    if ( params->d_db->keyExists( "Defaults" ) ) {
        auto defaults_db = params->d_db->getDatabase( "Defaults" ); // get defaults database
        auto defaultkeys = defaults_db->getAllKeys();               // get defaults database keys
        // check that all default argument names are in the list of property argument names map
        // each default string to a boolean indicating if the default was found in at least one
        // property argument list
        std::map<std::string, bool>
            defaults_found; // Defaults key to boolean for that key being found
        // initialize entries in defaults_found to false
        for ( auto &defaultkey : defaultkeys ) {
            defaults_found.insert( std::make_pair( defaultkey, false ) );
        }
        // for each property needed by formulation
        for ( auto &propertie : properties ) {
            // get its arguments
            std::vector<std::string> argnames =
                d_properties.find( propertie )->second->get_arguments();
            // for each Defaults key
            for ( auto &defaultkey : defaultkeys ) {
                // try to find it in the property arguments
                auto hit = std::find( argnames.begin(), argnames.end(), defaultkey );
                // if found, report it as being found
                if ( hit != argnames.end() )
                    defaults_found.find( defaultkey )->second = true;
            }
        }
        // generate error if a Defaults key was not found in any property arguments
        for ( auto &defaultkey : defaultkeys ) {
            std::string insist_string =
                "Default argument '" + defaultkey + "' was not found as a property argument";
            std::map<std::string, bool>::const_iterator it = defaults_found.find( defaultkey );
            bool found                                     = false;
            if ( it != defaults_found.end() ) {
                found = it->second;
            }
            AMP_INSIST( found, insist_string );
        }

        // load and check defaults:
        // for each property needed by formulation
        for ( auto &propertie : properties ) {
            auto property      = d_properties.find( propertie )->second; // pointer to property
            size_t n_arguments = property->get_number_arguments();       // number of arguments
            auto argnames      = property->get_arguments();              // argument names
            std::vector<double> prop_defaults( n_arguments );            // argument default values
            auto ranges = property->get_arg_ranges();                    // argument ranges
            // for each argument
            for ( size_t i = 0; i < n_arguments; ++i ) {
                // initially set default value to 1.0000001*(argument range minimum)
                prop_defaults[i] = ranges[i][0] * ( 1.0000001 );
                // try to find argument in Defaults keys
                auto hit = std::find( defaultkeys.begin(), defaultkeys.end(), argnames[i] );
                // if found,
                if ( hit != defaultkeys.end() ) {
                    // use the value provided with the Defaults key
                    prop_defaults[i] = defaults_db->getScalar<double>( argnames[i] );
                    // ensure that the default value is within the argument range
                    AMP_INSIST( property->in_range( argnames[i], prop_defaults[i] ),
                                std::string( "Default for argument " ) + argnames[i] +
                                    std::string( " is out of range" ) );
                } else {
                    AMP_WARNING( "Default value for key ''" + argnames[i] +
                                     "'' was not found in SubchannelPhysicsModel database. " +
                                     "Default value set to "
                                 << prop_defaults[i] );
                }
            }
            // set the defaults
            property->set_defaults( prop_defaults );
        }
    } else {
        AMP_ERROR( "Defaults database was not supplied in input database" );
    }
}

void SubchannelPhysicsModel::getProperty(
    std::string property,
    std::vector<double> &result,
    std::map<std::string, std::shared_ptr<std::vector<double>>> &args )
{
    // evaluate material property
    auto it = d_properties.find( property );
    AMP_INSIST( it != d_properties.end(), "Model does not have property (" + property + ")" );
    d_properties.find( property )->second->evalv( result, args );
}
} // namespace AMP::Operator
