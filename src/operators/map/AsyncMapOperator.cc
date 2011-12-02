#include "AsyncMapOperator.h"
#include "AsyncMapOperatorParameters.h"

namespace AMP {
namespace Operator {


AsyncMapOperator::AsyncMapOperator ( const boost::shared_ptr <OperatorParameters> &p )
    : AsynchronousOperator ( p )
{
    boost::shared_ptr<AsyncMapOperatorParameters> params = boost::dynamic_pointer_cast<AsyncMapOperatorParameters>(p);
    d_comm = params->d_MapComm;
    d_mesh1 = params->d_Mesh1;
    d_mesh2 = params->d_Mesh2;
}


AsyncMapOperator::~AsyncMapOperator ()
{
}


}
}

