#ifndef included_AMP_ParameterFactory
#define included_AMP_ParameterFactory

/* AMP files */
#include "AMP/mesh/Mesh.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/Database.h"

namespace AMP::Operator {

class ParameterFactory
{

public:
    ParameterFactory() {}
    ~ParameterFactory() {}

    static std::shared_ptr<OperatorParameters>
    createParameter( std::shared_ptr<AMP::Database> input_db,
                     std::shared_ptr<AMP::Mesh::Mesh> mesh );
};
} // namespace AMP::Operator

#endif
