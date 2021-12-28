#ifndef included_AMP_MassAssemblyParameters
#define included_AMP_MassAssemblyParameters

#include "AssemblyParameters.h"

namespace AMP::Operator {

class MassAssemblyParameters : public AssemblyParameters
{
public:
    explicit MassAssemblyParameters( std::shared_ptr<AMP::Database> db ) : AssemblyParameters( db )
    {
    }

    virtual ~MassAssemblyParameters() {}


protected:
private:
};
} // namespace AMP::Operator

#endif
