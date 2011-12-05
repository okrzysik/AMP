#include "AsyncMapOperator.h"
#include "AsyncMapColumnOperator.h"


namespace AMP {
namespace Operator {


size_t globalMapTagOffset = 0;      // Initialize the global map tag offset


AsyncMapColumnOperator::AsyncMapColumnOperator ( const boost::shared_ptr<OperatorParameters> & params )
    : AsynchronousColumnOperator ( params )
{
}


void  AsyncMapColumnOperator::setVector ( AMP::LinearAlgebra::Vector::shared_ptr &p )
{
    d_OutputVector = p;
    for (size_t i=0; i<d_Operators.size(); i++)
        boost::dynamic_pointer_cast<AsyncMapOperator>(d_Operators[i])->setVector ( d_OutputVector );
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


void AsyncMapColumnOperator::apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const  AMP::LinearAlgebra::Vector::shared_ptr &u, 
        AMP::LinearAlgebra::Vector::shared_ptr  &r,
        const double a, const double b)
{
    this->applyStart  ( f , u , r , a , b );
    this->applyFinish ( f , u , r , a , b );
    if ( requiresMakeConsistentSet() ) {
        AMP_ASSERT(d_OutputVector.get()!=NULL);
        d_OutputVector->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    }
}


bool AsyncMapColumnOperator::requiresMakeConsistentSet()
{ 
    bool test = false;
    for (size_t i=0; i<d_Operators.size(); i++)
        test = test | boost::dynamic_pointer_cast<AsyncMapOperator>(d_Operators[i])->requiresMakeConsistentSet();
    return test;
}


}
}

