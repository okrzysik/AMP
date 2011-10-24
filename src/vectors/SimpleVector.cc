#include "math.h"

#include "SimpleVector.h"


namespace AMP {
namespace LinearAlgebra {

  double SimpleVector::L1Norm(void) const
  {
    double ans = 0.0;
    for ( const_iterator cur = begin() ; cur != end() ; cur++ )
    {
      ans += fabs (*cur);
    }
    return ans;
  }

  double SimpleVector::L2Norm(void) const
  {
    double ans = 0.0;
    for ( const_iterator cur = begin() ; cur != end() ; cur++ )
    {
      ans += (*cur) * (*cur);
    }
    return sqrt ( ans );
  }

  double SimpleVector::maxNorm(void) const
  {
    double ans = 0.0;
    for ( const_iterator cur = begin() ; cur != end() ; cur++ )
    {
      ans = std::max ( ans , fabs (*cur) );
    }
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
    return ans;
  }

  void SimpleVector::setToScalar(double alpha)
  {
    for ( iterator cur = begin() ; cur != end() ; cur++ )
    {
      (*cur) = alpha;
    }
  }

  void SimpleVector::scale(double alpha, const VectorOperations &x)
  {
    iterator cur = begin();
    const_iterator curx = x.castTo<Vector>().begin();
    while ( cur != end() )
    {
      (*cur) = alpha * (*curx);
      cur++;
      curx++;
    }
  }

  void SimpleVector::setRandomValues(void)
  {
    srand ( time ( NULL ) );
    for ( iterator cur = begin() ; cur != end() ; cur++ )
    {
      (*cur) = (double)(rand() & 10000);
      if (*cur == 0.0 )
        *cur = 4.5;
      (*cur) /= 10000.0;
    }
  }

  void SimpleVector::scale(double alpha)
  {
    for ( iterator cur = begin() ; cur != end() ; cur++ )
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
      curx++; cury++; cur++;
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
      curx++; cury++; cur++;
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
      curx++; cury++; cur++;
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
      curx++; cury++; cur++;
    }
  }

  void SimpleVector::reciprocal(const VectorOperations &x)
  {
    iterator cur = begin();
    const_iterator curx = x.castTo<Vector>().begin();
    while ( cur != end() )
    {
      (*cur) = 1. / (*curx);
      cur++;
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
      curx++; cury++; cur++;
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
      curx++; cury++; cur++;
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
      curx++; cur++;
    }
  }

  void SimpleVector::abs(const VectorOperations &x)
  {
    iterator cur = begin();
    const_iterator curx = x.castTo<Vector>().begin();
    while ( cur != end() )
    {
      (*cur) = fabs (*curx);
      cur++;
      curx++;
    }
  }


}
}

