
#ifndef included_AMP_MechanicsModelParameters
#define included_AMP_MechanicsModelParameters

#include "AMP/operators/ElementPhysicsModelParameters.h"

#include "AMP/vectors/Vector.h"

namespace AMP {
namespace Operator {

/** A class encapsulating all the parameters that the material model
 * requires to evaluate the stress and/or tangent  */
class MechanicsModelParameters : public ElementPhysicsModelParameters
{
public:
    /** Constructor */
    explicit MechanicsModelParameters( AMP::shared_ptr<AMP::Database> db )
        : ElementPhysicsModelParameters( db )
    {
    }

    /** Destructor */
    virtual ~MechanicsModelParameters() {}

    /** A vector of deformation gradient values, which are required to
     * compute the stress and/or tangent. */
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_deformationGradient;
};
} // namespace Operator
} // namespace AMP

#endif
