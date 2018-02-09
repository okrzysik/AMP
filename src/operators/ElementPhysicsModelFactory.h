#ifndef included_AMP_ElementPhysicsModelFactory
#define included_AMP_ElementPhysicsModelFactory

/* Boost files */
#include "AMP/utils/shared_ptr.h"

#include "AMP/operators/ElementPhysicsModel.h"

namespace AMP {
namespace Operator {

class ElementPhysicsModelFactory
{
public:
    ElementPhysicsModelFactory() {}
    ~ElementPhysicsModelFactory() {}

    static AMP::shared_ptr<ElementPhysicsModel>
    createElementPhysicsModel( AMP::shared_ptr<AMP::Database> input_db );
};
} // namespace Operator
} // namespace AMP

#endif
