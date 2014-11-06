#ifndef included_AMP_Map3to1to3Parameters
#define included_AMP_Map3to1to3Parameters

#include "operators/map/AsyncMapOperatorParameters.h"


namespace AMP {
namespace Operator {


class Map3to1to3Parameters : public AsyncMapOperatorParameters
{
public:
    int          d_MasterValue;
    int          d_NumToSend;

    Map3to1to3Parameters ( const AMP::shared_ptr<AMP::Database> &db ):
        AsyncMapOperatorParameters(db), d_MasterValue(0), d_NumToSend(0)
    {
    }
};


}
}

#endif
