#include "math.h"

#include "SimpleVector.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Constructors                                                  *
****************************************************************/
SimpleVector::SimpleVector () : Vector ()
{
}
Vector::shared_ptr  SimpleVector::create ( size_t localSize , Variable::shared_ptr var )
{
    boost::shared_ptr<SimpleVector> retVal( new SimpleVector );
    retVal->d_startIndex = 0;
    retVal->setVariable ( var );
    retVal->d_Data.resize ( localSize );
    AMP_MPI comm(AMP_COMM_SELF);
    AMP::Discretization::DOFManager::shared_ptr DOFs( new AMP::Discretization::DOFManager( localSize, comm ) );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList( AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_comm = comm;
    retVal->d_globalSize = localSize;
    return retVal;
}
Vector::shared_ptr  SimpleVector::create ( size_t localSize , Variable::shared_ptr var, AMP_MPI comm )
{
    boost::shared_ptr<SimpleVector> retVal( new SimpleVector );
    retVal->d_startIndex = 0;
    retVal->setVariable ( var );
    retVal->d_Data.resize ( localSize );
    AMP::Discretization::DOFManager::shared_ptr DOFs( new AMP::Discretization::DOFManager( localSize, comm ) );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList( AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_comm = comm;
    retVal->d_globalSize = localSize;
    return retVal;
}
Vector::shared_ptr  SimpleVector::create ( Variable::shared_ptr var,
    AMP::Discretization::DOFManager::shared_ptr DOFs, 
    AMP::LinearAlgebra::CommunicationList::shared_ptr commlist )
{
    boost::shared_ptr<SimpleVector> retVal( new SimpleVector );
    retVal->d_startIndex = DOFs->beginDOF();
    retVal->setVariable ( var );
    retVal->d_Data.resize (  DOFs->numLocalDOF() );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList( commlist );
    retVal->d_comm = DOFs->getComm();
    retVal->d_globalSize = DOFs->numGlobalDOF();
    return retVal;
}

        

/****************************************************************
* Compute min, max, norms, and dot product                      *
****************************************************************/
double SimpleVector::min(void) const
{
    double local_min = *std::min_element ( d_Data.begin() , d_Data.end() );
    return d_comm.minReduce(local_min);
}
double SimpleVector::max(void) const
{
    double local_max =*std::max_element ( d_Data.begin() , d_Data.end() );
    return d_comm.maxReduce(local_max);
}
double SimpleVector::L1Norm(void) const
{
    double ans = 0.0;
    for ( const_iterator cur = begin() ; cur != end() ; ++cur )
        ans += fabs (*cur);
    ans = d_comm.sumReduce(ans);
    return ans;
}
double SimpleVector::L2Norm(void) const
{
    double ans = 0.0;
    for ( const_iterator cur = begin() ; cur != end() ; ++cur )
        ans += (*cur) * (*cur);
    ans = d_comm.sumReduce(ans);
    return sqrt ( ans );
}
double SimpleVector::maxNorm(void) const
{
    double ans = 0.0;
    for ( const_iterator cur = begin() ; cur != end() ; ++cur )
        ans = std::max ( ans , fabs (*cur) );
    ans = d_comm.maxReduce(ans);
    return ans;
}
double SimpleVector::dot(const VectorOperations &rhs ) const
{
    AMP_ASSERT ( getLocalSize() > 0 );
    AMP_ASSERT ( getLocalSize() == rhs.castTo<Vector>().getLocalSize() );
    const_iterator a = begin();
    const_iterator b = rhs.castTo<Vector>().begin();
    const_iterator a_end = end();
    //const_iterator b_end = rhs.castTo<Vector>().end();
    double  ans = 0;
    while ( a != a_end )
    {
      ans += *a * *b;
      b++;
      a++;
    }
    ans = d_comm.sumReduce(ans);
    return ans;
}


/****************************************************************
* Scale the vector and set to scalar                            *
****************************************************************/
void SimpleVector::setToScalar(double alpha)
{
    for ( iterator cur = begin() ; cur != end() ; ++cur )
        (*cur) = alpha;
    this->makeConsistent(CONSISTENT_SET);
}
void SimpleVector::scale(double alpha, const VectorOperations &x)
{
    iterator cur = begin();
    const_iterator curx = x.castTo<Vector>().begin();
    while ( cur != end() )
    {
      (*cur) = alpha * (*curx);
      ++cur;
      curx++;
    }
}
void SimpleVector::scale(double alpha)
{
    for ( iterator cur = begin() ; cur != end() ; ++cur )
    {
      (*cur) *= alpha;
    }
}



  void SimpleVector::add(const VectorOperations &x, const VectorOperations &y)
  {
    const_iterator curx , cury;
    iterator cur;
    curx = x.castTo<Vector>().begin();
    cury = y.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
      (*cur) = (*curx) + (*cury );
      curx++; cury++; ++cur;
    }
  }

  void SimpleVector::subtract(const VectorOperations &x, const VectorOperations &y)
  {
    const_iterator curx , cury; 
    iterator cur;
    curx = x.castTo<Vector>().begin();
    cury = y.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
      (*cur) = (*curx) - (*cury );
      curx++; cury++; ++cur;
    }
  }

  void SimpleVector::multiply( const VectorOperations &x, const VectorOperations &y)
  {
    const_iterator curx , cury;
    iterator cur;
    curx = x.castTo<Vector>().begin();
    cury = y.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
      (*cur) = (*curx) * (*cury );
      curx++; cury++; ++cur;
    }
  }

  void SimpleVector::divide( const VectorOperations &x, const VectorOperations &y)
  {
    const_iterator curx , cury;
    iterator cur;
    curx = x.castTo<Vector>().begin();
    cury = y.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
      (*cur) = (*curx) / (*cury );
      curx++; cury++; ++cur;
    }
  }

  void SimpleVector::reciprocal(const VectorOperations &x)
  {
    iterator cur = begin();
    const_iterator curx = x.castTo<Vector>().begin();
    while ( cur != end() )
    {
      (*cur) = 1. / (*curx);
      ++cur;
      curx++;
    }
  }

  void SimpleVector::linearSum(double alpha, const VectorOperations &x,
                      double beta, const VectorOperations &y)
  {
    const_iterator curx , cury;
    iterator cur;
    curx = x.castTo<Vector>().begin();
    cury = y.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
      (*cur) = alpha * (*curx) + beta * (*cury);
      curx++; cury++; ++cur;
    }
  }

  void SimpleVector::axpy(double alpha, const VectorOperations &x, const VectorOperations &y)
  {
    const_iterator curx , cury;
    iterator cur;
    curx = x.castTo<Vector>().begin();
    cury = y.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
      (*cur) = alpha * (*curx) + (*cury);
      curx++; cury++; ++cur;
    }
  }

  void SimpleVector::axpby(double alpha, double beta, const VectorOperations &x)
  {
    const_iterator curx;
    iterator cur;
    curx = x.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
      (*cur) = alpha * (*curx) + beta * (*cur);
      curx++; ++cur;
    }
  }

  void SimpleVector::abs(const VectorOperations &x)
  {
    iterator cur = begin();
    const_iterator curx = x.castTo<Vector>().begin();
    while ( cur != end() )
    {
      (*cur) = fabs (*curx);
      ++cur;
      curx++;
    }
  }


}
}

