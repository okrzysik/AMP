#ifndef included_ParameterFactory
#define included_ParameterFactory

/* AMP files */
#include "ampmesh/Mesh.h"
#include "operators/OperatorParameters.h"
#include "utils/Database.h"

namespace AMP {
namespace Operator {

class ParameterFactory
{

public:
    ParameterFactory() {}
    ~ParameterFactory() {}

    static AMP::shared_ptr<OperatorParameters>
    createParameter( AMP::shared_ptr<AMP::Database> input_db, AMP::Mesh::Mesh::shared_ptr mesh );
};
} // namespace Operator
} // namespace AMP

#endif
