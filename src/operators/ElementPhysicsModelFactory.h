#ifndef included_AMP_ElementPhysicsModelFactory
#define included_AMP_ElementPhysicsModelFactory

#include "AMP/operators/ElementPhysicsModel.h"
#include <memory>

namespace AMP::Operator {

class ElementPhysicsModelFactory
{
public:
    ElementPhysicsModelFactory() {}
    ~ElementPhysicsModelFactory() {}

    static std::shared_ptr<ElementPhysicsModel>
    createElementPhysicsModel( std::shared_ptr<AMP::Database> input_db );
};
} // namespace AMP::Operator

#endif
