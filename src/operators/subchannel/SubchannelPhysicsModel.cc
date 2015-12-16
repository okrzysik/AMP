/*
 * SubchannelPhysicsModel.cc
 *
 */

#include "operators/subchannel/SubchannelPhysicsModel.h"
#include "utils/Database.h"
#include "utils/Utilities.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>

namespace AMP {
namespace Operator {

SubchannelPhysicsModel::SubchannelPhysicsModel(
    const AMP::shared_ptr<ElementPhysicsModelParameters> &params )
    : ElementPhysicsModel( params )
{
    // get material key
    AMP_INSIST( ( params->d_db->keyExists( "Material" ) ),
                "Subchannel Key ''Material'' is missing!" );
    std::string matname = params->d_db->getString( "Material" );
    d_material = AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( matname );

    // get the formulation key
    AMP_INSIST( ( params->d_db->keyExists( "Formulation" ) ),
                "Subchannel Key ''Formulation'' is missing!" );
    std::string formulation = params->d_db->getString( "Formulation" );

    // determine which property functions are needed based on formulation key
    std::vector<std::string> properties;
    if ( formulation == std::string( "OneEqnForwardSubstitution" ) ) {
        AMP_ERROR( "The formulation ''OneEqnForwardSubstitution'' has not yet been implemented" );
    } else if ( formulation == std::string( "TwoEqnPicardIteration" ) ) {
        properties.push_back( "Enthalpy" );
    } else if ( formulation == std::string( "TwoEqnJFNKPressure" ) ) {
        properties.push_back( "Enthalpy" );
        properties.push_back( "SpecificVolume" );
    } else if ( formulation == std::string( "TwoEqnJFNKDensity" ) ) {
        properties.push_back( "Enthalpy" );
        properties.push_back( "Pressure" );
        AMP_ERROR( "The formulation ''TwoEqnJFNKDensity'' has not yet been implemented" );
    } else if ( formulation == std::string( "FunctionsTesting" ) ) {
        properties.push_back( "Temperature" );
        properties.push_back( "SaturatedLiquidEnthalpy" );
        properties.push_back( "SpecificVolume" );
        properties.push_back( "ThermalConductivity" );
        properties.push_back( "ConvectiveHeat" );
        properties.push_back( "DynamicViscosity" );
        properties.push_back( "Enthalpy" );
    } else {
        AMP_ERROR( "Invalid Formulation key" );
    }
    // add properties to property pointer map
    for ( std::vector<std::string>::iterator prop = properties.begin(); prop != properties.end();
          ++prop ) {
        d_properties.insert( std::make_pair( *prop, d_material->property( *prop ) ) );
    }

    if ( params->d_db->keyExists( "Defaults" ) ) {
        AMP::shared_ptr<Database> defaults_db =
            params->d_db->getDatabase( "Defaults" ); // get defaults database
        std::vector<std::string> defaultkeys =
            defaults_db->getAllKeys(); // get defaults database keys
        // check that all default argument names are in the list of property argument names
        // map each default string to a boolean indicating if the default was found in at least one
        // property argument
        // list
        std::map<std::string, bool>
            defaults_found; // maps Defaults key to boolean for that key being found
        // initialize entries in defaults_found to false
        for ( std::vector<std::string>::iterator key = defaultkeys.begin();
              key != defaultkeys.end();
              ++key ) {
            defaults_found.insert( std::make_pair( *key, false ) );
        }
        // for each property needed by formulation
        for ( std::vector<std::string>::iterator prop = properties.begin();
              prop != properties.end();
              ++prop ) {
            // get its arguments
            std::vector<std::string> argnames = d_properties.find( *prop )->second->get_arguments();
            // for each Defaults key
            for ( std::vector<std::string>::iterator key = defaultkeys.begin();
                  key != defaultkeys.end();
                  ++key ) {
                // try to find it in the property arguments
                std::vector<std::string>::iterator hit =
                    std::find( argnames.begin(), argnames.end(), *key );
                // if found, report it as being found
                if ( hit != argnames.end() )
                    defaults_found.find( *key )->second = true;
            }
        }
        // generate error if a Defaults key was not found in any property arguments
        for ( std::vector<std::string>::iterator key = defaultkeys.begin();
              key != defaultkeys.end();
              ++key ) {
            std::string insist_string =
                "Default argument '" + ( *key ) + "' was not found as a property argument";
            std::map<std::string, bool>::const_iterator it = defaults_found.find( *key );
            bool found = false;
            if ( it != defaults_found.end() ) {
                found = it->second;
            }
            AMP_INSIST( found, insist_string );
        }

        // load and check defaults:
        // for each property needed by formulation
        for ( std::vector<std::string>::iterator prop = properties.begin();
              prop != properties.end();
              ++prop ) {
            AMP::Materials::PropertyPtr property =
                d_properties.find( *prop )->second;                        // pointer to property
            size_t n_arguments = property->get_number_arguments();         // number of arguments
            std::vector<std::string> argnames = property->get_arguments(); // argument names
            std::vector<double> prop_defaults( n_arguments ); // argument default values
            std::vector<std::vector<double>> ranges = property->get_arg_ranges(); // argument ranges
            // for each argument
            for ( size_t i = 0; i < n_arguments; ++i ) {
                // initially set default value to 1.0000001*(argument range minimum)
                prop_defaults[i] = ranges[i][0] * ( 1.0000001 );
                // try to find argument in Defaults keys
                std::vector<std::string>::iterator hit =
                    std::find( defaultkeys.begin(), defaultkeys.end(), argnames[i] );
                // if found,
                if ( hit != defaultkeys.end() ) {
                    // use the value provided with the Defaults key
                    prop_defaults[i] = defaults_db->getDouble( argnames[i] );
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
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> &args )
{
    // evaluate material property
    std::map<std::string, AMP::Materials::PropertyPtr>::iterator it = d_properties.find( property );
    AMP_INSIST( it != d_properties.end(), "Model does not have property (" + property + ")" );
    d_properties.find( property )->second->evalv( result, args );
}
}
}
