#include "AMP/operators/diffusion/DiffusionTransportTensorModel.h"
#include "AMP/utils/Database.h"

#include <cmath>


namespace AMP::Operator {

DiffusionTransportTensorModel::DiffusionTransportTensorModel(
    std::shared_ptr<const DiffusionTransportTensorModelParameters> params )
    : DiffusionTransportModel( params )
{
    d_IsTensor = true;

    std::string modelName = params->d_db->getString( "name" );
    if ( not( modelName == "DiffusionCylindricalTransportModel" ) ) {
        AMP_INSIST( d_property->isTensor(), "material property must be of tensor type" );
    }
}

void DiffusionTransportTensorModel::getTensorTransport(
    AMP::Array<std::shared_ptr<std::vector<double>>> &result,
    std::map<std::string, std::shared_ptr<std::vector<double>>> &args,
    const std::vector<libMesh::Point> & )
{
    std::shared_ptr<std::vector<double>> scaledp;
    double lower, upper;

    if ( d_UseBilogScaling ) {
        // do the transform
        lower                       = d_BilogRange[0] + d_BilogEpsilonRangeLimit;
        upper                       = d_BilogRange[1] - d_BilogEpsilonRangeLimit;
        scaledp                     = bilogTransform( ( *args[d_BilogVariable] ), lower, upper );
        std::vector<double> &scaled = *scaledp;

        // save untransformed argument value
        for ( size_t i = 0; i < ( *args[d_BilogVariable] ).size(); i++ ) {
            double temp                   = ( *args[d_BilogVariable] )[i];
            ( *args[d_BilogVariable] )[i] = scaled[i];
            scaled[i]                     = temp;
        }
    }

    // evaluate material property
    // material library has been temporarily supplied with a dummy evalv for tensors
    // new material interface will fix.
    d_property->evalv( result, {}, args );

    if ( d_UseBilogScaling ) {
        // restore untransformed argument value
        std::vector<double> &scaled = *scaledp;
        lower                       = d_BilogRange[0] + d_BilogEpsilonRangeLimit;
        upper                       = d_BilogRange[1] - d_BilogEpsilonRangeLimit;
        for ( size_t i = 0; i < ( *args[d_BilogVariable] ).size(); i++ ) {
            ( *args[d_BilogVariable] )[i] = scaled[i];
        }

        if ( d_BilogScaleCoefficient ) {
            for ( auto &elem : result )
                bilogScale( *elem, lower, upper );
        }
    }
}
} // namespace AMP::Operator
