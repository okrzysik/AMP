#include "math.h"

#include "SimpleVector.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Constructors                                                  *
****************************************************************/
template <typename T>
SimpleVector<T>::SimpleVector (): 
    Vector(),
    d_startIndex(0),
    d_globalSize(0)
{
}

template <typename T>
Vector::shared_ptr  SimpleVector<T>::create ( size_t localSize , Variable::shared_ptr var )
{
    AMP::shared_ptr<SimpleVector<T> > retVal( new SimpleVector<T> );
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

template <typename T>
Vector::shared_ptr  SimpleVector<T>::create ( size_t localSize , Variable::shared_ptr var, AMP_MPI comm )
{
    AMP::shared_ptr<SimpleVector<T> > retVal( new SimpleVector<T> );
    retVal->d_startIndex = 0;
    retVal->setVariable ( var );
    retVal->d_Data.resize ( localSize );
    AMP::Discretization::DOFManager::shared_ptr DOFs( new AMP::Discretization::DOFManager( localSize, comm ) );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList( AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_comm = comm;
    retVal->d_globalSize = comm.sumReduce(localSize);
    return retVal;
}

template <typename T>
Vector::shared_ptr  SimpleVector<T>::create ( Variable::shared_ptr var,
    AMP::Discretization::DOFManager::shared_ptr DOFs, 
    AMP::LinearAlgebra::CommunicationList::shared_ptr commlist )
{
    AMP::shared_ptr<SimpleVector<T> > retVal( new SimpleVector<T> );
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
template <typename T>
double SimpleVector<T>::min(void) const
{
    auto local_min = static_cast<double> (*std::min_element ( d_Data.begin() , d_Data.end() ) );
    return d_comm.minReduce(local_min);
}

template <typename T>
double SimpleVector<T>::max(void) const
{
    auto local_max =static_cast<double> (*std::max_element ( d_Data.begin() , d_Data.end() ));
    return d_comm.maxReduce(local_max);
}

template <typename T>
double SimpleVector<T>::L1Norm(void) const
{
    double ans = 0.0;

    for ( const T &val:d_Data )
        ans += static_cast<double> (fabs (val));

    ans = d_comm.sumReduce(ans);
    return ans;
}

template <typename T>
double SimpleVector<T>::L2Norm(void) const
{
    double ans = 0.0;
    for ( const auto &val:d_Data )
        ans += static_cast<double> (val * val);
    ans = d_comm.sumReduce(ans);
    return sqrt ( ans );
}

template <typename T>
double SimpleVector<T>::maxNorm(void) const
{
    double ans = 0.0;
    for ( const auto &val:d_Data )
        ans = static_cast<double> (std::max ( ans , fabs (val) ));
    ans = d_comm.maxReduce(ans);
    return ans;
}

template <typename T>
double SimpleVector<T>::dot(const VectorOperations &rhs ) const
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
        ans += static_cast<double> (*a * *b);
      ++b;
      ++a;
    }
    ans = d_comm.sumReduce(ans);
    return ans;
}


/****************************************************************
* Copy vector                                                   *
****************************************************************/
template <typename T>
void SimpleVector<T>::copyVector( Vector::const_shared_ptr src_vec )
{
    if ( getLocalSize() != src_vec->getLocalSize() )
        AMP_ERROR( "Mismatched vectors" );
    ConstVectorDataIterator it = src_vec->begin();
    for (size_t i=0; i<getLocalSize(); i++) {
        d_Data[i] = static_cast<T> (*it);
        ++it;
    }
    copyGhostValues( src_vec );
    // Copy the consistency state from the rhs
    *d_UpdateState = *(src_vec->getUpdateStatusPtr());
}


/****************************************************************
* Copy raw data                                                 *
****************************************************************/
template <typename T>
void SimpleVector<T>::putRawData ( const double *in )
{
    for (size_t i=0; i<d_Data.size(); ++i) {
        d_Data[i] = static_cast<T>(in[i]);
    }
}

template <typename T>
void SimpleVector<T>::copyOutRawData ( double *out ) const
{
    for (size_t i=0; i<d_Data.size(); ++i) {
        out[i] = static_cast<double>(d_Data[i]);
    }
}


/****************************************************************
* Scale the vector and set to scalar                            *
****************************************************************/
template <typename T>
void SimpleVector<T>::setToScalar(double alpha)
{
    auto alpha_T = static_cast<T>(alpha);

    for ( auto &val: d_Data )
        val = alpha_T;
    this->makeConsistent(CONSISTENT_SET);
}

template <typename T>
void SimpleVector<T>::scale(double alpha, const VectorOperations &x)
{
    auto alpha_T = static_cast<T>(alpha);

    iterator cur = begin();
    const_iterator curx = x.castTo<Vector>().begin();
    while ( cur != end() )
    {
        (*cur) = alpha_T * static_cast<T>((*curx));
        ++cur;
        ++curx;
    }
}

template <typename T>
void SimpleVector<T>::scale(double alpha)
{
    auto alpha_T = static_cast<T>(alpha);

    for ( auto &val:d_Data )
    {
        val *= alpha_T;
    }
}

template <typename T>
void SimpleVector<T>::add(const VectorOperations &x, const VectorOperations &y)
{
    const_iterator curx , cury;
    iterator cur;
    curx = x.castTo<Vector>().begin();
    cury = y.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
        (*cur) = static_cast<T>(*curx) + static_cast<T>(*cury );
        ++curx; ++cury; ++cur;
    }
}

template <typename T>
void SimpleVector<T>::subtract(const VectorOperations &x, const VectorOperations &y)
{
    const_iterator curx , cury; 
    iterator cur;
    curx = x.castTo<Vector>().begin();
    cury = y.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
        (*cur) = static_cast<T>(*curx) - static_cast<T>(*cury );
        ++curx; ++cury; ++cur;
    }
}

template <typename T>
void SimpleVector<T>::multiply( const VectorOperations &x, const VectorOperations &y)
{
    const_iterator curx , cury;
    iterator cur;
    curx = x.castTo<Vector>().begin();
    cury = y.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
        (*cur) = static_cast<T>(*curx) * static_cast<T>(*cury );
        ++curx; ++cury; ++cur;
    }
}

template <typename T>
void SimpleVector<T>::divide( const VectorOperations &x, const VectorOperations &y)
{
    const_iterator curx , cury;
    iterator cur;
    curx = x.castTo<Vector>().begin();
    cury = y.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
        (*cur) = static_cast<T>(*curx) / static_cast<T>(*cury );
        ++curx; ++cury; ++cur;
    }
}

template <typename T>
void SimpleVector<T>::reciprocal(const VectorOperations &x)
{
    iterator cur = begin();
    const_iterator curx = x.castTo<Vector>().begin();
    while ( cur != end() )
    {
        (*cur) = 1. / static_cast<T>(*curx);
        ++cur;
        ++curx;
    }
}

template <typename T>
void SimpleVector<T>::linearSum(double alpha, const VectorOperations &x,
                  double beta, const VectorOperations &y)
{
    const_iterator curx , cury;
    iterator cur;
    curx = x.castTo<Vector>().begin();
    cury = y.castTo<Vector>().begin();
    cur = begin();
    auto alpha_T  = static_cast<T>(alpha);
    auto beta_T  = static_cast<T>(beta);
    while ( cur != end() )
    {
        (*cur) = alpha_T * static_cast<T>(*curx) + beta_T * static_cast<T>(*cury);
        ++curx; ++cury; ++cur;
    }
}

template <typename T>
void SimpleVector<T>::axpy(double alpha, const VectorOperations &x, const VectorOperations &y)
{
    const_iterator curx , cury;
    iterator cur;
    auto alpha_T  = static_cast<T>(alpha);
    curx = x.castTo<Vector>().begin();
    cury = y.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
        (*cur) = alpha_T * static_cast<T>(*curx) + static_cast<T>(*cury);
        ++curx; ++cury; ++cur;
    }
}

template <typename T>
void SimpleVector<T>::axpby(double alpha, double beta, const VectorOperations &x)
{
    const_iterator curx;
    iterator cur;
    auto alpha_T  = static_cast<T>(alpha);
    auto beta_T  = static_cast<T>(beta);
    curx = x.castTo<Vector>().begin();
    cur = begin();
    while ( cur != end() )
    {
        (*cur) = alpha_T * static_cast<T>(*curx) + beta_T * static_cast<T>(*cur);
        ++curx; ++cur;
    }
}

template <typename T>
void SimpleVector<T>::abs(const VectorOperations &x)
{
    iterator cur = begin();
    const_iterator curx = x.castTo<Vector>().begin();
    while ( cur != end() )
    {
        (*cur) = fabs(static_cast<T>(*curx));
        ++cur;
        ++curx;
    }
}


}
}

