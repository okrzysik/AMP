#ifndef included_AMP_FlowTransportModel
#define included_AMP_FlowTransportModel

#include <cstring>

#include "AMP/materials/Material.h"
#include "AMP/operators/ElementPhysicsModel.h"
#include <memory>


namespace AMP {
namespace Operator {


typedef ElementPhysicsModelParameters FlowTransportModelParameters;

class FlowTransportModel : public ElementPhysicsModel
{
public:
    explicit FlowTransportModel( std::shared_ptr<FlowTransportModelParameters> params )
        : ElementPhysicsModel( params )
    {
        d_useMaterialsLibrary =
            params->d_db->getWithDefault<bool>( "USE_MATERIALS_LIBRARY", false );

        if ( d_useMaterialsLibrary == true ) {
            AMP_INSIST( ( params->d_db->keyExists( "Material" ) ), "Key ''Material'' is missing!" );
            std::string matname = params->d_db->getString( "Material" );
            d_coolant =
                AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( matname );
        } else {
            d_density = params->d_db->getWithDefault<double>( "DENSITY", 1 );
            d_fmu     = params->d_db->getWithDefault<double>( "VISCOSITY", 1.0 );
            d_Re      = params->d_db->getWithDefault<double>( "ReynoldsNumber", 1.0 );
        }
    }


    /**
     * Destructor.
     */
    virtual ~FlowTransportModel() {}

    double getDensity() { return d_density; }

    double getViscosity() { return d_fmu; }

    double getReynoldsNumber() { return d_Re; }

protected:
    bool d_useMaterialsLibrary;

    double d_density, d_fmu, d_Re;

    std::shared_ptr<AMP::Materials::Material>
        d_coolant; /**< Shared pointer to the materials object. */

private:
};
} // namespace Operator
} // namespace AMP

#endif
