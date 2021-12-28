
#ifndef included_AMP_BlockOperatorParameters
#define included_AMP_BlockOperatorParameters

#include "AMP/operators/OperatorParameters.h"

#include <memory>

namespace AMP::Operator {

class BlockOperatorParameters : public OperatorParameters
{
public:
    explicit BlockOperatorParameters( std::shared_ptr<AMP::Database> db ) : OperatorParameters( db )
    {
    }

    virtual ~BlockOperatorParameters() {}

    std::vector<std::vector<std::shared_ptr<OperatorParameters>>> d_blockParams;
};
} // namespace AMP::Operator

#endif
