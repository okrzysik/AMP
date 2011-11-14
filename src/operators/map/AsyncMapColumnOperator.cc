#include "AsyncMapOperator.h"
#include "AsyncMapColumnOperator.h"


namespace AMP {
namespace Operator {


  AsyncMapColumnOperator::AsyncMapColumnOperator ( const boost::shared_ptr<OperatorParameters> & params )
    : AsynchronousColumnOperator ( params )
  {
  }

  void  AsyncMapColumnOperator::setVector ( AMP::LinearAlgebra::Vector::shared_ptr &p )
  {
    std::vector< boost::shared_ptr < Operator > >::iterator curOp = d_Operators.begin();
    while ( curOp != d_Operators.end() )
    {
      boost::dynamic_pointer_cast<AsyncMapOperator> ( *curOp )->setVector ( p );
      curOp++;
    }
  }

  void  AsyncMapColumnOperator::append ( boost::shared_ptr < Operator > op )
  {
    boost::shared_ptr<AsyncMapColumnOperator>  mapColumn = boost::dynamic_pointer_cast<AsyncMapColumnOperator> ( op );
    if ( mapColumn )
    {
      std::vector< boost::shared_ptr < Operator > >::iterator curOp = mapColumn.get()->d_Operators.begin();
      while ( curOp != mapColumn.get()->d_Operators.end() )
      {
        append ( *curOp );
        curOp++;
      }
    }
    else
    {
      boost::shared_ptr<AsyncMapOperator>  mapOp = boost::dynamic_pointer_cast<AsyncMapOperator> ( op );
      AMP_INSIST ( mapOp , "Attempt to add a non-AsyncMapOperator to a AsyncMapColumnOperator" );
      AsynchronousColumnOperator::append ( mapOp );
    }
  }

}
}

