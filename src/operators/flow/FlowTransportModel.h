
#ifndef included_AMP_FlowTransportModel
#define included_AMP_FlowTransportModel

#include <cstring>

#include "operators/ElementPhysicsModel.h"

#include "materials/Material.h"

#include "utils/shared_ptr.h"

namespace AMP {
namespace Operator {

typedef ElementPhysicsModelParameters FlowTransportModelParameters;

class FlowTransportModel : public ElementPhysicsModel
{
public:
    explicit FlowTransportModel( const AMP::shared_ptr<FlowTransportModelParameters> &params )
        : ElementPhysicsModel( params )
    {
        d_useMaterialsLibrary =
            ( params->d_db )->getBoolWithDefault( "USE_MATERIALS_LIBRARY", false );

        if ( d_useMaterialsLibrary == true ) {
            AMP_INSIST( ( params->d_db->keyExists( "Material" ) ), "Key ''Material'' is missing!" );
            std::string matname = params->d_db->getString( "Material" );
            d_coolant =
                AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( matname );
        } else {
            d_density = ( params->d_db )->getDoubleWithDefault( "DENSITY", 1 );
            d_fmu     = ( params->d_db )->getDoubleWithDefault( "VISCOSITY", 1.0 );
            d_Re      = ( params->d_db )->getDoubleWithDefault( "ReynoldsNumber", 1.0 );
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

    AMP::Materials::Material::shared_ptr d_coolant; /**< Shared pointer to the materials object. */

private:
};
} // namespace Operator
} // namespace AMP

#endif
