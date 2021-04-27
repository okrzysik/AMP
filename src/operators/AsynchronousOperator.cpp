

#include "AsynchronousOperator.h"

namespace AMP {
namespace Operator {


AsynchronousOperator::AsynchronousOperator( std::shared_ptr<const OperatorParameters> params )
    : Operator( params )
{
    // Initialize the request list to an empty vector
    d_RequestList.clear();
}


AsynchronousOperator::~AsynchronousOperator() = default;


/*bool AsynchronousOperator::continueAsynchronousConstruction( std::shared_ptr<const
OperatorParameters>& )
{
    return true;
}*/


void AsynchronousOperator::waitForAllRequests()
{
    auto curReq = beginRequests();
    while ( curReq != endRequests() ) {
        MPI_Request request = *curReq;
        if ( request != 0 )
            AMP_MPI::wait( request );
        ++curReq;
    }
}


void AsynchronousOperator::clearRequests() { d_RequestList.clear(); }


void AsynchronousOperator::reserveRequests( size_t i ) { d_RequestList.resize( i ); }


std::vector<MPI_Request>::iterator AsynchronousOperator::beginRequests()
{
    return d_RequestList.begin();
}


std::vector<MPI_Request>::iterator AsynchronousOperator::endRequests()
{
    return d_RequestList.end();
}


MPI_Request &AsynchronousOperator::getRequest( size_t i )
{
    AMP_ASSERT( i < d_RequestList.size() );
    return d_RequestList[i];
}

void AsynchronousOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                  AMP::LinearAlgebra::Vector::shared_ptr f )
{
    applyStart( u, f );
    applyFinish( u, f );
}
} // namespace Operator
} // namespace AMP
