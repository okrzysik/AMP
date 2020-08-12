#ifndef included_AMP_Vector_inline
#define included_AMP_Vector_inline

#include "AMP/vectors/data/VectorDataIterator.h"
#include <algorithm>


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Get basic info                                                *
 ****************************************************************/
inline std::shared_ptr<ParameterBase> Vector::getParameters()
{
    return std::shared_ptr<ParameterBase>();
}
inline AMP::Discretization::DOFManager::shared_ptr Vector::getDOFManager() const
{
    return d_DOFManager;
}
inline AMP_MPI Vector::getComm() const { return d_CommList->getComm(); }


/****************************************************************
 * Subset for variable name                                      *
 ****************************************************************/
inline Vector::shared_ptr Vector::subsetVectorForVariable( const std::string &name )
{
    auto var = std::make_shared<AMP::LinearAlgebra::Variable>( name );
    return subsetVectorForVariable( var );
}
inline Vector::const_shared_ptr
Vector::constSubsetVectorForVariable( const std::string &name ) const
{
    auto var = std::make_shared<AMP::LinearAlgebra::Variable>( name );
    return constSubsetVectorForVariable( var );
}


/****************************************************************
 * getView/hasView                                               *
 ****************************************************************/
template<typename VIEW_TYPE>
Vector::shared_ptr Vector::getView() const
{
    for ( size_t i = 0; i != d_Views->size(); i++ ) {
        auto vec = ( *d_Views )[i].lock();
        if ( vec ) {
            auto vec2 = std::dynamic_pointer_cast<VIEW_TYPE>( vec );
            if ( vec2 ) {
                return Vector::shared_ptr( ( *d_Views )[i] );
            }
        }
    }
    return Vector::shared_ptr();
}
template<typename VIEW_TYPE>
bool Vector::hasView() const
{
    for ( size_t i = 0; i != d_Views->size(); i++ ) {
        auto vec = ( *d_Views )[i].lock();
        if ( vec ) {
            auto vec2 = std::dynamic_pointer_cast<VIEW_TYPE>( vec );
            if ( vec2 ) {
                return true;
            }
        }
    }
    return false;
}


/****************************************************************
 * RNG                                                           *
 ****************************************************************/
inline void Vector::setDefaultRNG( RNG::shared_ptr p ) { d_DefaultRNG = p; }
inline RNG::shared_ptr Vector::getDefaultRNG()
{
    if ( !d_DefaultRNG ) {
        AMP_MPI globalComm( AMP_COMM_WORLD );
        int rank = globalComm.getRank();
        RNGParameters::shared_ptr params( new RNGParameters(
            RNGParameters::RNGOptions::USE_GLOBAL_SEED, static_cast<size_t>( rank ) ) );
        d_DefaultRNG = RNG::shared_ptr( new RNG( params ) );
    }
    return d_DefaultRNG;
}


/****************************************************************
 * Misc functions                                                *
 ****************************************************************/
inline const Variable::shared_ptr Vector::getVariable() const { return d_pVariable; }
inline Variable::shared_ptr Vector::getVariable() { return d_pVariable; }
inline Vector::shared_ptr Vector::cloneVector() const { return cloneVector( getVariable() ); }
inline void Vector::setVariable( const Variable::shared_ptr name )
{
    AMP_ASSERT( name.get() != nullptr );
    d_pVariable = name;
}


/****************************************************************
 * Wrappers for shared_ptr                                       *
 ****************************************************************/
// clang-format off
inline void Vector::swapVectors( shared_ptr other ) { swapVectors( *other ); }
inline void Vector::aliasVector( shared_ptr other ) { aliasVector( *other ); }
inline void Vector::addCommunicationListToParameters( CommunicationList::shared_ptr ) {}
inline void Vector::aliasGhostBuffer( shared_ptr in ) { d_Ghosts = in->d_Ghosts; }
inline std::ostream &operator<<( std::ostream &out, const Vector::shared_ptr p ) { return operator<<( out, *p ); }
// clang-format on

/****************************************************************
 * Math API for Vector                                          *
 ****************************************************************/

inline void Vector::copy( const VectorData &x )
{
  d_VectorOps->copy(x, *d_VectorData);
}

inline void Vector::zero( void )
{
  d_VectorOps->zero(*d_VectorData);
}

inline void Vector::setToScalar( double alpha )
{
  d_VectorOps->setToScalar(alpha, *d_VectorData);
}

inline void Vector::setRandomValues( void )
{
  d_VectorOps->setRandomValues(*d_VectorData);
}

inline void Vector::setRandomValues( RNG::shared_ptr rng )
{
  d_VectorOps->setRandomValues(rng, *d_VectorData);
}

inline void Vector::scale( double alpha, const VectorData &x )
{
  d_VectorOps->scale(alpha, x, *d_VectorData);
}

inline void Vector::scale( double alpha )
{
  d_VectorOps->scale(alpha, *d_VectorData);
}

inline void Vector::add( const VectorData &x, const VectorData &y )
{
  d_VectorOps->add(x, y, *d_VectorData);
}

inline void Vector::subtract( const VectorData &x, const VectorData &y )
{
  d_VectorOps->subtract(x, y, *d_VectorData);
}

inline void Vector::multiply( const VectorData &x, const VectorData &y )
{
  d_VectorOps->multiply(x, y, *d_VectorData);
}

inline void Vector::divide( const VectorData &x, const VectorData &y )
{
  d_VectorOps->divide(x, y, *d_VectorData);
}

inline void Vector::reciprocal( const VectorData &x )
{
  d_VectorOps->reciprocal(x, *d_VectorData);
}

inline void Vector::linearSum( double alpha, const VectorData &x, double beta, const VectorData &y )
{
  d_VectorOps->linearSum(alpha, x, beta, y, *d_VectorData);
}

inline void Vector::axpy( double alpha, const VectorData &x, const VectorData &y )
{
  d_VectorOps->axpy(alpha, x, y, *d_VectorData);
}

inline void Vector::axpby( double alpha, double beta, const VectorData &x )
{
  d_VectorOps->axpby(alpha, beta, x, *d_VectorData);
}

inline void Vector::abs( const VectorData &x )
{
  d_VectorOps->abs(x, *d_VectorData);
}

inline void Vector::addScalar( const VectorData &x, double alpha_in )
{
  d_VectorOps->addScalar(x, alpha_in, *d_VectorData);
}

inline double Vector::min( void ) const
{
  return d_VectorOps->min(*d_VectorData);
}

inline double Vector::max( void ) const
{
  return d_VectorOps->max(*d_VectorData);
}

inline double Vector::L1Norm( void ) const
{
  std::cout << "Entering Vector::L1Norm " << std::endl;
  std::cout << "Calling L1Norm on " << typeid(d_VectorData).name() << " with VecOps " << typeid(d_VectorOps).name() << std::endl;
  return d_VectorOps->L1Norm(*d_VectorData);
  std::cout << "Exiting Vector::L1Norm " << std::endl;
}

inline double Vector::L2Norm( void ) const
{
  return d_VectorOps->L2Norm(*d_VectorData);
}

inline double Vector::maxNorm( void ) const
{
  return d_VectorOps->maxNorm(*d_VectorData);
}

inline double Vector::minQuotient( const VectorData &x ) const
{
  return d_VectorOps->minQuotient(x, *d_VectorData);
}

inline double Vector::wrmsNorm( const VectorData &x, const VectorData &y ) const
{
  return d_VectorOps->wrmsNorm(x, y);
}

inline double Vector::wrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y ) const
{
  return d_VectorOps->wrmsNormMask(x, mask, y);
}

inline double Vector::dot( const VectorData &x ) const
{
  return d_VectorOps->dot(x, *d_VectorData);
}
    
inline bool Vector::equals( const VectorData &a, double tol ) const
{
  return d_VectorOps->equals(a, *d_VectorData, tol);
}

inline double Vector::localMin( void ) const
{
  return d_VectorOps->localMin(*d_VectorData);
}

inline double Vector::localMax( void ) const
{
  return d_VectorOps->localMax(*d_VectorData);
}

inline double Vector::localL1Norm( void ) const
{
  return d_VectorOps->localL1Norm(*d_VectorData);
}

inline double Vector::localL2Norm( void ) const
{
  return d_VectorOps->localL2Norm(*d_VectorData);
}

inline double Vector::localMaxNorm( void ) const
{
  return d_VectorOps->localMaxNorm(*d_VectorData);
}

inline double Vector::localDot( const VectorData &x ) const
{
  return d_VectorOps->localDot(x, *d_VectorData);
}

inline double Vector::localMinQuotient( const VectorData &x ) const
{
  return d_VectorOps->localMinQuotient(x, *d_VectorData);
}

inline double Vector::localWrmsNorm( const VectorData &x ) const
{
  return d_VectorOps->localWrmsNorm(x, *d_VectorData);
}

inline double Vector::localWrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y ) const
{
  return d_VectorOps->localWrmsNormMask(x, mask, y);
}

inline bool Vector::localEquals( const VectorData &x, double tol ) const
{
  return d_VectorOps->localEquals(x, *d_VectorData, tol);
}

/****************************************************************
 * Wrappers for shared_ptr                                       *
 ****************************************************************/

inline bool Vector::equals( std::shared_ptr<const VectorData> x,
			    double tol ) const
{
  return equals( *x, tol );
}

inline void Vector::copy( std::shared_ptr<const VectorData> x )
{
    copy( *x );
}

inline void Vector::scale( double alpha,
			   std::shared_ptr<const VectorData> x )
{
    scale( alpha, *x );
}

inline void Vector::add( std::shared_ptr<const VectorData> x,
			 std::shared_ptr<const VectorData> y )
{
    add( *x, *y );
}

void Vector::addScalar( std::shared_ptr<const VectorData> x,
			double alpha )
{
    addScalar( *x, alpha );
}

inline void Vector::subtract( std::shared_ptr<const VectorData> x,
			      std::shared_ptr<const VectorData> y )
{
    subtract( *x, *y );
}
inline void Vector::multiply( std::shared_ptr<const VectorData> x,
			      std::shared_ptr<const VectorData> y )
{
    multiply( *x, *y );
}
inline void Vector::divide( std::shared_ptr<const VectorData> x,
			    std::shared_ptr<const VectorData> y )
{
    divide( *x, *y );
}
inline void Vector::reciprocal( std::shared_ptr<const VectorData> x )
{
    reciprocal( *x );
}
inline void Vector::linearSum( double alpha,
			       std::shared_ptr<const VectorData> x,
			       double beta,
			       std::shared_ptr<const VectorData> y )
{
    linearSum( alpha, *x, beta, *y );
}
inline void Vector::axpy( double alpha,
			  std::shared_ptr<const VectorData> x,
			  std::shared_ptr<const VectorData> y)
{
    axpy( alpha, *x, *y );
}
inline void Vector::axpby( double alpha,
			   double beta,
			   std::shared_ptr<const VectorData> x )
{
    axpby( alpha, beta, *x );
}
inline void Vector::abs( std::shared_ptr<const VectorData> x )
{
    return abs( *x );
}
inline double Vector::dot( std::shared_ptr<const VectorData> x ) const
{
    return dot( *x );
}

inline double Vector::minQuotient( std::shared_ptr<const VectorData> x ) const
{
    return minQuotient( *x );
}
inline double Vector::wrmsNorm( std::shared_ptr<const VectorData> x,
				std::shared_ptr<const VectorData> y ) const
{
    return d_VectorOps->wrmsNorm( *x, *y );
}
inline double Vector::wrmsNormMask( std::shared_ptr<const VectorData> x,
				    std::shared_ptr<const VectorData> mask,
				    std::shared_ptr<const VectorData> y ) const
{
    return d_VectorOps->wrmsNormMask( *x, *mask, *y );
}
   
} // namespace LinearAlgebra
} // namespace AMP

#endif
