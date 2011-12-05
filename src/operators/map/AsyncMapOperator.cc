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


void AsyncMapOperator::apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const  AMP::LinearAlgebra::Vector::shared_ptr &u, 
        AMP::LinearAlgebra::Vector::shared_ptr  &r,
        const double a, const double b)
{
    applyStart  ( f , u , r , a , b );
    applyFinish ( f , u , r , a , b );
    if ( requiresMakeConsistentSet() ) {
        AMP_ASSERT(d_OutputVector.get()!=NULL);
        d_OutputVector->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }
}


bool AsyncMapOperator::requiresMakeConsistentSet()
{ 
    return false;
}


}
}

