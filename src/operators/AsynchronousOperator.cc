

#include "AsynchronousOperator.h"

namespace AMP {
namespace Operator {


AsynchronousOperator::AsynchronousOperator
                       ( const AMP::shared_ptr < OperatorParameters > &params )
      : Operator ( params )
{
    // Initialize the request list to an empty vector
    d_RequestList.clear ();
}


AsynchronousOperator::~AsynchronousOperator ()
{
}


/*bool AsynchronousOperator::continueAsynchronousConstruction ( const AMP::shared_ptr < OperatorParameters > & )
{
    return true;
}*/


void AsynchronousOperator::waitForAllRequests ()
{
    std::vector<MPI_Request>::iterator  curReq = beginRequests ();
    while ( curReq != endRequests() )
    {
        MPI_Request request = *curReq;
        if ( request!=0 )
            AMP_MPI::wait( request );
        ++curReq;
    }
}


void AsynchronousOperator::clearRequests ()
{
    d_RequestList.clear ();
}


void AsynchronousOperator::reserveRequests ( size_t i )
{
    d_RequestList.resize ( i );
}


std::vector<MPI_Request>::iterator  AsynchronousOperator::beginRequests ()
{
    return d_RequestList.begin();
}


std::vector<MPI_Request>::iterator  AsynchronousOperator::endRequests ()
{
    return d_RequestList.end();
}


MPI_Request &AsynchronousOperator::getRequest ( size_t i )
{
    AMP_ASSERT ( i < d_RequestList.size() );
    return d_RequestList[i];
}

void AsynchronousOperator::apply(AMP::LinearAlgebra::Vector::const_shared_ptr f,
         AMP::LinearAlgebra::Vector::const_shared_ptr u, AMP::LinearAlgebra::Vector::shared_ptr r,
         const double a , const double b )
{
    applyStart ( f , u , r , a , b );
    applyFinish ( f , u , r , a , b );
}



}
}
