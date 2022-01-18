#ifndef included_AMP_ElementOperationFactory
#define included_AMP_ElementOperationFactory

#include "AMP/operators/ElementOperation.h"
#include <memory>


namespace AMP::Operator {

class ElementOperationFactory
{
public:
    ElementOperationFactory() {}
    ~ElementOperationFactory() {}

    static std::shared_ptr<ElementOperation>
    createElementOperation( std::shared_ptr<AMP::Database> input_db );
};
} // namespace AMP::Operator

#endif
