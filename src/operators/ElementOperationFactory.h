#ifndef included_AMP_ElementOperationFactory
#define included_AMP_ElementOperationFactory

#include "AMP/operators/ElementOperation.h"
#include "AMP/utils/shared_ptr.h"


namespace AMP {
namespace Operator {

class ElementOperationFactory
{
public:
    ElementOperationFactory() {}
    ~ElementOperationFactory() {}

    static AMP::shared_ptr<ElementOperation>
    createElementOperation( AMP::shared_ptr<AMP::Database> input_db );
};
} // namespace Operator
} // namespace AMP

#endif
