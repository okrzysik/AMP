#ifndef included_AMP_MassAssemblyParameters
#define included_AMP_MassAssemblyParameters

#include "AssemblyParameters.h"

namespace AMP {
namespace Operator {

class MassAssemblyParameters : public AssemblyParameters
{
public:
    explicit MassAssemblyParameters( AMP::shared_ptr<AMP::Database> db ) : AssemblyParameters( db )
    {
    }

    virtual ~MassAssemblyParameters() {}


protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif
