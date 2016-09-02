/*
 * DiffusionCylindricalTransportModel.cc
 *
 *  Created on: Aug 19, 2011
 *      Author: gad
 */

#include "DiffusionCylindricalTransportModel.h"
#include "utils/Database.h"
#include "utils/Utilities.h"
#include <cmath>
#include <map>

/* Libmesh files */
#include "libmesh/point.h"

namespace AMP {
namespace Operator {

DiffusionCylindricalTransportModel::DiffusionCylindricalTransportModel(
    const AMP::shared_ptr<DiffusionCylindricalTransportModelParameters> params )
    : DiffusionTransportTensorModel( params )
{
    AMP_INSIST( params->d_db->keyExists( "RadiusArgument" ),
                "must specify RadiusArgument for material" );
    d_RadiusArgument = params->d_db->getString( "RadiusArgument" );
}

void DiffusionCylindricalTransportModel::getTensorTransport(
    std::vector<std::vector<AMP::shared_ptr<std::vector<double>>>> &result,
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> &args,
    const std::vector<libMesh::Point> &coordinates )
{
    AMP::shared_ptr<std::vector<double>> scaledp;
    double lower, upper;

    if ( d_UseBilogScaling ) {
        // do the transform
        lower                       = d_BilogRange[0] + d_BilogEpsilonRangeLimit;
        upper                       = d_BilogRange[1] - d_BilogEpsilonRangeLimit;
        scaledp                     = bilogTransform( *args[d_BilogVariable], lower, upper );
        std::vector<double> &scaled = *scaledp;

        // save untransformed argument value
        for ( size_t i = 0; i < args[d_BilogVariable]->size(); i++ ) {
            double temp                   = ( *args[d_BilogVariable] )[i];
            ( *args[d_BilogVariable] )[i] = scaled[i];
            scaled[i]                     = temp;
        }
    }

    // evaluate material property as a function of radius
    // first fill in radius array
    std::vector<double> radialCoefficient( coordinates.size() );
    if ( args.find( d_RadiusArgument ) != args.end() ) {
        AMP_INSIST( args[d_RadiusArgument]->size() == coordinates.size(), "radial size mismatch" );
    } else {
        args.insert( std::make_pair( d_RadiusArgument,
                                     AMP::shared_ptr<std::vector<double>>(
                                         new std::vector<double>( coordinates.size() ) ) ) );
    }
    std::vector<double> &radius( *args[d_RadiusArgument] );
    for ( size_t k = 0; k < radius.size(); k++ ) {
        double x  = coordinates[k]( 0 );
        double y  = coordinates[k]( 1 );
        double r  = sqrt( x * x + y * y );
        radius[k] = r;
    }
    d_property->evalv( radialCoefficient, args );

    if ( d_UseBilogScaling ) {
        // restore untransformed argument value
        std::vector<double> &scaled = *scaledp;
        lower                       = d_BilogRange[0] + d_BilogEpsilonRangeLimit;
        upper                       = d_BilogRange[1] - d_BilogEpsilonRangeLimit;
        for ( size_t i = 0; i < ( *args[d_BilogVariable] ).size(); i++ ) {
            ( *args[d_BilogVariable] )[i] = scaled[i];
        }

        if ( d_BilogScaleCoefficient ) {
            for ( size_t i = 0; i < radialCoefficient.size(); i++ ) {
                bilogScale( radialCoefficient, lower, upper );
            }
        }
    }

    // form angle-dependent tensor factor
    for ( size_t i = 0; i < 3; i++ ) {
        AMP_INSIST( result.size() == 3, "result tensor must be 3x3" );
        for ( size_t j = 0; j < 3; j++ ) {
            AMP_INSIST( result[i].size() == 3, "result tensor must be 3x3" );
            AMP_INSIST( result[i][j]->size() == coordinates.size(),
                        "result tensor components must be same size as coordinates" );
            if ( ( i == 2 ) || ( j == 2 ) )
                for ( size_t k           = 0; k < coordinates.size(); k++ )
                    ( *result[i][j] )[k] = 0.;
        }
    }
    for ( size_t k = 0; k < coordinates.size(); k++ ) {
        double x             = coordinates[k]( 0 );
        double y             = coordinates[k]( 1 );
        double r2            = x * x + y * y;
        ( *result[0][0] )[k] = radialCoefficient[k] * x * x / r2;
        ( *result[0][1] )[k] = radialCoefficient[k] * x * y / r2;
        ( *result[1][0] )[k] = radialCoefficient[k] * x * y / r2;
        ( *result[1][1] )[k] = radialCoefficient[k] * y * y / r2;
    }
}
}
}
