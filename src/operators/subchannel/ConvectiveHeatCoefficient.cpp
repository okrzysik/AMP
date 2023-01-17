#include "AMP/operators/subchannel/ConvectiveHeatCoefficient.h"

#include <cmath>


namespace AMP::Operator {


ConvectiveHeatCoefficient::ConvectiveHeatCoefficient(
    std::shared_ptr<const RobinPhysicsModelParameters> params )
    : RobinPhysicsModel( params )
{
    AMP_INSIST( params->d_db->keyExists( "Material" ),
                "Convective Heat  Coefficient Key ''Material'' is missing!" );
    auto matname = params->d_db->getString( "Material" );

    d_material = AMP::Materials::getMaterial( matname );

    AMP_INSIST( params->d_db->keyExists( "Property" ),
                "Convective Heat Coefficient Key ''Property'' is missing!" );
    auto propname = params->d_db->getString( "Property" );
    d_property    = d_material->property( propname );

    d_defaults = d_property->get_defaults();
    if ( params->d_db->keyExists( "Defaults" ) ) {
        // check for correct names
        auto defaults_db = params->d_db->getDatabase( "Defaults" );
        auto defaultkeys = defaults_db->getAllKeys();
        AMP_INSIST( defaultkeys.size() == d_property->get_number_arguments(),
                    "Incorrect number of defaults supplied." );
        d_argNames = d_property->get_arguments();
        for ( auto &defaultkey : defaultkeys ) {
            auto hit = std::find( d_argNames.begin(), d_argNames.end(), defaultkey );
            AMP_INSIST( hit != d_argNames.end(),
                        std::string( "Argument name " ) + defaultkey +
                            std::string( " is invalid" ) );
        }

        // load defaults into the material property, checking range validity
        for ( size_t i = 0; i < d_argNames.size(); ++i ) {
            d_defaults[i] = defaults_db->getScalar<double>( d_argNames[i] );
            AMP_INSIST( d_property->in_range( d_argNames[i], d_defaults[i] ),
                        std::string( "Default for argument " ) + d_argNames[i] +
                            std::string( " is out of range" ) );
        }
    }
    d_property->set_defaults( d_defaults );
}


void ConvectiveHeatCoefficient::getConductance(
    std::vector<double> &beta,
    std::vector<double> &gamma,
    const std::vector<std::vector<double>> &inputVectors )
{


    std::map<std::string, const std::vector<double> &> args = { { "temperature", inputVectors[0] },
                                                                { "density", inputVectors[2] },
                                                                { "diameter", inputVectors[3] } };
    //  { "reynolds", inputVectors[4] }, { "prandtl" , inputVectors[5] } };

    d_property->evalv( beta, {}, args );
    d_property->evalv( gamma, {}, args );
}
} // namespace AMP::Operator
