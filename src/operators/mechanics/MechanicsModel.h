#ifndef included_AMP_MechanicsModel
#define included_AMP_MechanicsModel

#include <cstring>

#include "AMP/materials/Material.h"
#include "AMP/operators/ElementPhysicsModel.h"
#include "AMP/operators/MechanicsModelParameters.h"
#include "AMP/vectors/Vector.h"
#include <memory>


namespace AMP {
namespace Operator {

class MechanicsModel : public ElementPhysicsModel
{
public:
    /** Constructor */
    explicit MechanicsModel( std::shared_ptr<const MechanicsModelParameters> &params )
        : ElementPhysicsModel( params )
    {
        bool useMaterialsLibrary = params->d_db->getWithDefault( "USE_MATERIALS_LIBRARY", false );
        if ( useMaterialsLibrary == true ) {
            AMP_INSIST( ( params->d_db->keyExists( "Material" ) ), "Key ''Material'' is missing!" );
            std::string matname = params->d_db->getString( "Material" );
            d_material =
                AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( matname );
        }
    }

    /** Destructor */
    virtual ~MechanicsModel() {}

    /** This function is used by the Mechanics operator to pass relevant
     * values to this material model. These values will be used to compute
     * the stress and/or tangent. */
    virtual void reset( std::shared_ptr<MechanicsModelParameters> params )
    {
        d_deformationGradient = params->d_deformationGradient;
    }

    /** This function should return the consistent/continuum tangent values in the vector
     * provided.   */
    virtual void getTangent( std::shared_ptr<AMP::LinearAlgebra::Vector> tangent ) {}

    /** This function should return the stress values in the vector
     * provided.   */
    virtual void getStress( std::shared_ptr<AMP::LinearAlgebra::Vector> stress ) {}

protected:
    /** The material object that can be used to get material properties */
    std::shared_ptr<AMP::Materials::Material> d_material;

    /** Pointer to the deformation gradient that is required to compute the
     * stress and/or tangent. */
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_deformationGradient;
};
} // namespace Operator
} // namespace AMP

#endif
