#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/materials/CylindricallySymmetric.h"
#include "AMP/materials/ScalarProperty.h"
#include "AMP/operators/diffusion/DiffusionTransportTensorModel.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <cmath>
#include <map>


template<class TYPE>
static inline bool is( std::shared_ptr<AMP::Materials::Property> prop )
{
    return std::dynamic_pointer_cast<TYPE>( prop ).get() != nullptr;
}


namespace AMP::Operator {

const std::vector<libMesh::Point> DiffusionTransportModel::d_DummyCoords =
    std::vector<libMesh::Point>( 0 );

DiffusionTransportModel::DiffusionTransportModel(
    std::shared_ptr<const DiffusionTransportModelParameters> params )
    : ElementPhysicsModel( params ), d_defaults( Diffusion::NUMBER_VARIABLES ), d_IsTensor( false )
{
    AMP_INSIST( params->d_db->keyExists( "Material" ), "Diffusion Key ''Material'' is missing!" );
    std::string matname = params->d_db->getString( "Material" );

    d_material = AMP::Materials::getMaterial( matname );

    AMP_INSIST( params->d_db->keyExists( "Property" ), "Diffusion Key ''Property'' is missing!" );
    std::string propname = params->d_db->getString( "Property" );
    d_property           = d_material->property( propname );

    // load and check defaults
    // initially set them to the minimum of the range plus a bit
    d_defaults.resize( d_property->get_number_arguments() );
    auto ranges = d_property->get_arg_ranges();
    for ( size_t i = 0; i < d_defaults.size(); ++i ) {
        d_defaults[i] = ranges[i][0] * ( 1.0000001 );
    }
    if ( params->d_db->keyExists( "Defaults" ) ) {
        // check for correct names
        auto defaults_db = params->d_db->getDatabase( "Defaults" );
        auto defaultkeys = defaults_db->getAllKeys();
        AMP_INSIST( defaultkeys.size() == d_property->get_number_arguments(),
                    "Incorrect number of defaults supplied." );
        auto argnames = d_property->get_arguments();
        for ( auto &defaultkey : defaultkeys ) {
            auto hit = std::find( argnames.begin(), argnames.end(), defaultkey );
            AMP_INSIST( hit != argnames.end(),
                        std::string( "Argument name " ) + defaultkey +
                            std::string( " is invalid" ) );
        }

        // load defaults into the material property, checking range validity
        for ( size_t i = 0; i < argnames.size(); ++i ) {
            d_defaults[i] = defaults_db->getScalar<double>( argnames[i] );
            AMP_INSIST( d_property->in_range( argnames[i], d_defaults[i] ),
                        std::string( "Default for argument " ) + argnames[i] +
                            std::string( " is out of range" ) );
        }
    }
    d_property->set_defaults( d_defaults );

    // process bilog scaling details
    d_UseBilogScaling = params->d_db->getWithDefault<bool>( "UseBilogScaling", false );
    if ( d_UseBilogScaling ) {
        AMP_INSIST( params->d_db->keyExists( "BilogVariable" ), "must specify BilogVariable" );
        d_BilogVariable = params->d_db->getWithDefault<std::string>( "BilogVariable", "NONE" );

        d_BilogRange = d_property->get_arg_range( d_BilogVariable );
        AMP_INSIST( d_BilogRange[1] > d_BilogRange[0],
                    "material argument upper bound <= lower bound" );

        d_BilogScaleCoefficient =
            params->d_db->getWithDefault<bool>( "BilogScaleCoefficient", true );
        d_BilogEpsilonRangeLimit =
            params->d_db->getWithDefault<double>( "BilogEpsilonRangeLimit", 1.0e-06 );
    }

    // set property parameters
    if ( params->d_db->keyExists( "Parameters" ) ) {
        auto name     = d_property->get_name();
        auto data     = params->d_db->getVector<double>( "Parameters" );
        auto units    = d_property->get_units();
        auto args     = d_property->get_arguments();
        auto ranges   = d_property->get_arg_ranges();
        auto argUnits = d_property->get_arg_units();
        if ( is<AMP::Materials::CylindricallySymmetricTensor>( d_property ) ) {
            d_property =
                std::make_shared<AMP::Materials::CylindricallySymmetricTensor>( name, data );
        } else if ( is<AMP::Materials::PolynomialProperty>( d_property ) ) {
            d_property = std::make_shared<AMP::Materials::PolynomialProperty>(
                name, "", units, data, args, ranges, argUnits );
        } else if ( propname == "ThermalDiffusionCoefficient" && data.size() == 2 ) {
            d_property =
                std::make_shared<AMP::Materials::ScalarProperty>( name, data[0] * data[1], units );
        } else if ( d_property->isTensor() ) {
            auto dims = d_property->size();
            if ( params->d_db->keyExists( "Dimensions" ) )
                dims = params->d_db->getVector<size_t>( "Dimensions" );
            AMP_INSIST( dims.size() == 2, "only two dimensions allowed for tensor property" );
            AMP_ASSERT( data.size() == dims[0] * dims[1] );
            d_property =
                std::make_shared<AMP::Materials::ScalarTensorProperty>( name, "", dims, data );
        } else if ( d_property->isVector() ) {
            AMP_ASSERT( data.size() == 1 );
            d_property = std::make_shared<AMP::Materials::ScalarVectorProperty>( name, data[0] );
        } else {
            d_property = std::make_shared<AMP::Materials::PolynomialProperty>(
                name, "", units, data, args, ranges, argUnits );
        }
    }
}

std::shared_ptr<std::vector<double>> DiffusionTransportModel::bilogTransform(
    const std::vector<double> &U, const double a, const double b )
{
    auto up                = std::make_shared<std::vector<double>>( U.size() );
    std::vector<double> &u = *up;

    for ( size_t i = 0; i < U.size(); i++ ) {
        double eU = std::exp( U[i] );
        u[i]      = ( a + b * eU ) / ( 1.0 + eU );
    }

    return up;
}

void DiffusionTransportModel::bilogScale( std::vector<double> &v, const double a, const double b )
{
    for ( auto &elem : v ) {
        double ev     = std::exp( elem );
        double temp   = ( 1.0 + ev );
        double factor = ( b - a ) * ev / ( temp * temp );
        elem *= factor;
    }
}

void DiffusionTransportModel::getTransport(
    std::vector<double> &result,
    std::map<std::string, std::shared_ptr<std::vector<double>>> &args,
    const std::vector<libMesh::Point> & )
{
    PROFILE_START( "getTransport", 7 );
    std::shared_ptr<std::vector<double>> scaledp;


    if ( d_UseBilogScaling ) {
        // do the transform
        auto &data   = *args[d_BilogVariable];
        auto lower   = d_BilogRange[0] + d_BilogEpsilonRangeLimit;
        auto upper   = d_BilogRange[1] - d_BilogEpsilonRangeLimit;
        scaledp      = bilogTransform( data, lower, upper );
        auto &scaled = *scaledp;

        // save untransformed argument value
        for ( size_t i = 0; i < data.size(); i++ ) {
            double temp = data[i];
            data[i]     = scaled[i];
            scaled[i]   = temp;
        }
    }

    // evaluate material property
    d_property->evalv( result, args );

    if ( d_UseBilogScaling ) {
        // restore untransformed argument value
        auto &data                  = *args[d_BilogVariable];
        std::vector<double> &scaled = *scaledp;
        auto lower                  = d_BilogRange[0] + d_BilogEpsilonRangeLimit;
        auto upper                  = d_BilogRange[1] - d_BilogEpsilonRangeLimit;
        for ( size_t i = 0; i < data.size(); i++ ) {
            data[i] = scaled[i];
        }

        if ( d_BilogScaleCoefficient )
            bilogScale( result, lower, upper );
    }
    PROFILE_STOP( "getTransport", 7 );
}


} // namespace AMP::Operator
