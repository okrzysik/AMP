#include "math.h"

#include "ArrayVector.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Constructors                                                  *
****************************************************************/
template <typename T, typename FUN, typename Allocator>
ArrayVector<T, FUN, Allocator>::ArrayVector() : Vector()
{
}

template <typename T, typename FUN, typename Allocator>
Vector::shared_ptr ArrayVector<T, FUN, Allocator>::create( const std::vector<size_t> &localSize,
                                           Variable::shared_ptr var )
{
    AMP::shared_ptr<ArrayVector<T,FUN,Allocator>> retVal( new ArrayVector<T,FUN,Allocator>() );
    retVal->setVariable( var );

    retVal->resize( localSize );
    const auto N = retVal->getArray().length();
    AMP_MPI comm( AMP_COMM_SELF );
    AMP::Discretization::DOFManager::shared_ptr DOFs(
        new AMP::Discretization::DOFManager( N, comm ) );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_globalSize = N;
    retVal->d_comm       = comm;

    return retVal;
}

template <typename T, typename FUN, typename Allocator>
Vector::shared_ptr ArrayVector<T, FUN, Allocator>::create( const std::vector<size_t> &localSize,
                                           Variable::shared_ptr var,
                                           AMP_MPI comm )
{
    AMP::shared_ptr<ArrayVector<T,FUN,Allocator>> retVal( new ArrayVector<T,FUN,Allocator>() );
    retVal->setVariable( var );
    retVal->resize( localSize );
    const auto N = retVal->getArray().length();

    AMP::Discretization::DOFManager::shared_ptr DOFs(
        new AMP::Discretization::DOFManager( N, comm ) );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_comm       = comm;
    retVal->d_globalSize = comm.sumReduce( N );
    return retVal;
}

template <typename T, typename FUN, typename Allocator>
Vector::shared_ptr
ArrayVector<T, FUN, Allocator>::create( Variable::shared_ptr var,
                        AMP::Discretization::DOFManager::shared_ptr DOFs,
                        AMP::LinearAlgebra::CommunicationList::shared_ptr commlist )
{
    AMP::shared_ptr<ArrayVector<T,FUN,Allocator>> retVal( new ArrayVector<T,FUN,Allocator>() );
    retVal->setVariable( var );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList( commlist );
    retVal->d_comm       = DOFs->getComm();
    AMP_ERROR("This routine is not complete");
    return retVal;
}

/****************************************************************
* Copy vector                                                   *
****************************************************************/
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::copyVector( Vector::const_shared_ptr src_vec )
{
    if ( getLocalSize() != src_vec->getLocalSize() )
        AMP_ERROR( "Mismatched vectors" );
    auto &internalArray = this->getArray();
    const auto &otherArray = std::dynamic_pointer_cast<const ArrayVector<T, FUN, Allocator>>(src_vec)->getArray();
    internalArray = otherArray;  // use copy assignment
    // copy ghost values will also copy the consistency state
    copyGhostValues( src_vec );
}

template <typename T, typename FUN, typename Allocator>
inline Vector::shared_ptr ArrayVector<T, FUN, Allocator>::cloneVector( const Variable::shared_ptr name ) const
{
    const auto &array = this->getArray();
    return create( array.size(), name, this->getComm() );
}

template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::swapVectors( Vector &rhs )
{
    // get internal arrays
    AMP::Array<T, FUN, Allocator> &internalArray = this->getArray();
    AMP::Array<T, FUN, Allocator> &otherArray = dynamic_cast<ArrayVector<T, FUN, Allocator> &>(rhs).getArray();
    // reset views
    internalArray.swap( otherArray );
}

template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::aliasVector( Vector & )
{
    AMP_ERROR( "Not implemented" );
}

template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::resize( const std::vector<size_t> &localDims )
{
    d_array.resize(localDims);
}

// overloads for the VectorOperations that should go away eventually
    
/**
 * \param  alpha a scalar double
 * \brief  Set all compenents of a vector to a scalar.
 * For Vectors, the components of <em>this</em> are set to \f$\alpha\f$.
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::setToScalar( double alpha )
{
    T alpha_t = static_cast<T>(alpha);
    auto &array = this->getArray();
    array.fill(alpha_t);
}
/**
 * \param  alpha  a scalar double
 * \param  x  a vector
 * \brief  Set vector equal to scaled input.
 * For Vectors, \f$\mathit{this}_i = \alpha x_i\f$.
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::scale( double alpha, const VectorOperations &x )
{
    T alpha_t = static_cast<T>(alpha);
    const auto &xa = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(x);
    const auto &xarray = xa.getArray();
    auto fun = [&alpha_t](const T &xval){ return alpha_t*xval; };
    auto &array = this->getArray();
    FUN::transform(fun, xarray, array);
}

/**
 * \param  alpha  a scalar double
 *
 * \brief  Scale a vector.
 * For Vectors, \f$\mathit{this}_i = \alpha\mathit{this}_i\f$.
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::scale( double alpha )
{
    auto alpha_t = static_cast<T>(alpha);
    auto &array = this->getArray();
    array.scale(alpha_t);
}

/**
 * \param  x  a vector
 * \param  y  a vector
 * \brief  Adds two vectors.
 * For Vectors, \f$\mathit{this}_i = x_i + y_i\f$.
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::add( const VectorOperations &x, const VectorOperations &y )
{
    const auto &xa = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(x);
    const auto &ya = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(y);
    auto fun = []( const T &xval, const T &yval){ return xval+yval; };
    auto &array = this->getArray();
    const auto &xarray = xa.getArray();
    const auto &yarray = ya.getArray();
    FUN::transform(fun, xarray, yarray, array);
}

/**
 * \param x  a vector
 * \param y  a vector
 * \brief Subtracts one vector from another.
 * For Vectors, \f$\mathit{this}_i = x_i - y_i\f$
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::subtract( const VectorOperations &x, const VectorOperations &y )
{
    const auto &xa = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(x);
    const auto &ya = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(y);
    auto fun = []( const T &xval, const T &yval){ return xval-yval; };
    auto &array = this->getArray();
    const auto &xarray = xa.getArray();
    const auto &yarray = ya.getArray();
    FUN::transform(fun, xarray, yarray, array);
}

/**
 * \param x  a vector
 * \param y  a vector
 * \brief Component-wise multiply one vector with another.
 * For Vectors, \f$\mathit{this}_i = x_i  y_i\f$
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::multiply( const VectorOperations &x, const VectorOperations &y )
{
    const auto &xa = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(x);
    const auto &ya = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(y);
    auto fun = []( const T &xval, const T &yval){ return xval*yval; };
    auto &array = this->getArray();
    const auto &xarray = xa.getArray();
    const auto &yarray = ya.getArray();
    FUN::transform(fun, xarray, yarray, array);
}

/**
 * \param x  a vector
 * \param y  a vector
 * \brief Component-wise divide one vector by another.
 * For Vectors, \f$\mathit{this}_i = x_i / y_i\f$
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::divide( const VectorOperations &x, const VectorOperations &y )
{
    const auto &xa = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(x);
    const auto &ya = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(y);
    auto fun = []( const T &xval, const T &yval){ return xval/yval; };
    auto &array = this->getArray();
    const auto &xarray = xa.getArray();
    const auto &yarray = ya.getArray();
    FUN::transform(fun, xarray, yarray, array);
}

/**
 * \param x  a vector
 * \brief Set this to the component-wise reciprocal of a vector.  \f$\mathit{this}_i =
 * 1/x_i\f$.
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::reciprocal( const VectorOperations &x )
{
    const auto &xa = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(x);
    auto fun = []( const T &xval){ return 1.0/xval; };
    auto &array = this->getArray();
    const auto &xarray = xa.getArray();
    FUN::transform(fun, xarray, array);
}


/**
 * \param alpha a scalar
 * \param x a vector
 * \param beta a scalar
 * \param y a vector
 * \brief Set a vector to be a linear combination of two vectors.
 *  \f$\mathit{this}_i = \alpha x_i + \beta y_i\f$.
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::linearSum( double alpha,
                                                 const VectorOperations &x,
                                                 double beta,
                                                 const VectorOperations &y )
{
    auto alpha_t = static_cast<T>(alpha);
    auto beta_t = static_cast<T>(beta);
    const auto &xa = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(x);
    const auto &ya = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(y);
    auto fun = [&alpha_t, &beta_t]( const T &xval, const T &yval){ return alpha_t*xval+beta_t*yval; };
    auto &array = this->getArray();
    const auto &xarray = xa.getArray();
    const auto &yarray = ya.getArray();
    FUN::transform(fun, xarray, yarray, array);
}

/**
 * \param alpha a scalar
 * \param x a vector
 * \param y a vector
 * \brief Set this vector to alpha * x + y.  \f$\mathit{this}_i = \alpha x_i + y_i\f$.
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::axpy( double alpha, const VectorOperations &x, const VectorOperations &y )
{
    auto alpha_t = static_cast<T>(alpha);
    const auto &xa = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(x);
    const auto &ya = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(y);
    auto fun = [&alpha_t]( const T &xval, const T &yval){ return alpha_t*xval+yval; };
    auto &array = this->getArray();
    const auto &xarray = xa.getArray();
    const auto &yarray = ya.getArray();
    FUN::transform(fun, xarray, yarray, array);
}

/**
 * \param alpha a scalar
 * \param beta a scalar
 * \param x  a vector
 * \brief Set this vector alpha * x + this.
 * \f$\mathit{this}_i = \alpha x_i + \beta \mathit{this}_i \f$
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::axpby( double alpha, double beta, const VectorOperations &x )
{
    const auto alpha_t = static_cast<T>(alpha);
    const auto beta_t = static_cast<T>(beta);
    const auto &xa = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(x);
    const auto &xarray = xa.getArray();
    AMP::Array<T, FUN, Allocator> &array = this->getArray();
    array.axpby(alpha_t, xarray, beta_t);
}

/**
 * \param x a vector
 * \brief Set this to the component-wise absolute value of a vector.
 * \f$\mathit{this}_i = |x_i|\f$.
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::abs( const VectorOperations &x )
{
    const auto &xa = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(x);
    auto fun = []( const T &xval){ return std::abs(xval); };
    auto &array = this->getArray();
    const auto &xarray = xa.getArray();
    FUN::transform(fun, xarray, array);
}

/**
 * \brief Return the minimum value of the vector.  \f$\min_i \mathit{this}_i\f$.
 */
template <typename T, typename FUN, typename Allocator>
double ArrayVector<T, FUN, Allocator>::min( void ) const
{
    auto &array = this->getArray();
    auto local_min = array.min();
    auto global_min = d_comm.minReduce( local_min );
    return global_min;
}

/**
 * \brief Return the maximum value of the vector.  \f$\max_i \mathit{this}_i\f$.
 */
template <typename T, typename FUN, typename Allocator>
double ArrayVector<T, FUN, Allocator>::max( void ) const
{
    auto &array = this->getArray();
    auto local_max = array.max();
    auto global_max = d_comm.maxReduce( local_max );
    return global_max;
}

/**
 * \brief Set data in this vector to random values on [0,1).
 */
template <typename T, typename FUN, typename Allocator>
void ArrayVector<T, FUN, Allocator>::setRandomValues( void )
{
    auto &array = this->getArray();
    array.rand();
}

/**
 * \brief Return discrete @f$ L_1 @f$ -norm of this vector.
 * \details Returns \f[\sum_i |\mathit{this}_i|\f]
 */
template <typename T, typename FUN, typename Allocator>
double ArrayVector<T, FUN, Allocator>::L1Norm( void ) const
{
    auto fun = []( const T &xval, T &accum){ return accum+std::abs(xval); };
    auto &array = this->getArray();
    auto local_norm = FUN::sum(fun, array);
    auto global_norm = d_comm.sumReduce( local_norm );
    return global_norm;
}

/**
 * \brief Return discrete @f$ L_2 @f$ -norm of this vector.
 * \details Returns \f[\sqrt{\sum_i \mathit{this}_i^2}\f]
 */
template <typename T, typename FUN, typename Allocator>
double ArrayVector<T, FUN, Allocator>::L2Norm( void ) const
{
    auto fun = []( const T &xval, T &accum){ return accum+xval*xval; };
    auto &array = this->getArray();
    auto local_norm = FUN::sum(fun, array);
    auto global_norm = d_comm.sumReduce( local_norm );
    return std::sqrt(global_norm);
}

/**
 * \brief Return the @f$ L_\infty @f$ -norm of this vector.
 * \details Returns \f[\max_i |\mathit{this}_i|\f]
 */
template <typename T, typename FUN, typename Allocator>
double ArrayVector<T, FUN, Allocator>::maxNorm( void ) const
{
    auto fun = []( const T &xval, T &tmax){ return std::max(std::abs(tmax),std::abs(xval)); };
    auto &array = this->getArray();
    auto local_norm = FUN::reduce(fun, array);
    auto global_norm = d_comm.maxReduce( local_norm );
    return global_norm;
}

/**
 * \param x a vector
 * \brief Return the dot product of this vector with the argument vector.
 * \details Returns \f[\sum_i x_i\mathit{this}_i\f]
 */
template <typename T, typename FUN, typename Allocator>
double ArrayVector<T, FUN, Allocator>::dot( const VectorOperations &x ) const
{
    const auto &xa = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(x);
    auto fun = []( const T &xval, const T &yval, T &dotp){ return dotp+xval*yval; };
    auto &array = this->getArray();
    const auto &xarray = xa.getArray();
    auto local_dotp = FUN::reduce(fun, xarray, array);
    auto global_dotp = d_comm.sumReduce( local_dotp );
    return global_dotp;
}

template <typename T, typename FUN, typename Allocator>
bool ArrayVector<T, FUN, Allocator>::equals( Vector const &x, double tol ) const
{
    const auto &xa = dynamic_cast<const ArrayVector<T, FUN, Allocator> &>(x);
    const auto &xarray = xa.getArray();
    const auto &array = this->getArray();

    auto local_eq = array.equals(xarray, tol);
    auto global_eq = d_comm.allReduce(local_eq);
    return global_eq;
}

template <typename T, typename FUN, typename Allocator>
void  ArrayVector<T, FUN, Allocator>::putRawData( const double *buf )
{
    auto &array = this->getArray();
    array.copy(buf);
}

template <typename T, typename FUN, typename Allocator>
void  ArrayVector<T, FUN, Allocator>::copyOutRawData( double *buf ) const
{
    auto &array = this->getArray();
    array.copyTo(buf);
}

}
}
