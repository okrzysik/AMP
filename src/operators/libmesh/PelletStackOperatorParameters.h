
#ifndef included_AMP_PelletStackOperatorParameters
#define included_AMP_PelletStackOperatorParameters

#include "AMP/operators/map/AsyncMapColumnOperator.h"

namespace AMP {
namespace Operator {

class PelletStackOperatorParameters : public OperatorParameters
{
public:
    explicit PelletStackOperatorParameters( AMP::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
        d_currentPellet = static_cast<unsigned int>( -1 );
    }

    virtual ~PelletStackOperatorParameters() {}

    unsigned int d_currentPellet;
    AMP_MPI d_pelletStackComm;
    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator> d_n2nMaps;
};
} // namespace Operator
} // namespace AMP


#endif
