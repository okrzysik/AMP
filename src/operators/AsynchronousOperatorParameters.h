#ifndef included_AMP_AsynchronousOperatorParameters
#define included_AMP_AsynchronousOperatorParameters

#include "AMP/operators/OperatorParameters.h"

namespace AMP::Operator {


class AsynchronousOperatorParameters : public OperatorParameters
{
public:
    char d_AsynchronousConstructionParam;
    size_t d_ConstructionPhase;

    explicit AsynchronousOperatorParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db ), d_AsynchronousConstructionParam( 0 ), d_ConstructionPhase( 0 )
    {
    }
};
} // namespace AMP::Operator

#endif
