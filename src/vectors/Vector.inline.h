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
inline void Vector::aliasGhostBuffer( shared_ptr in ) { d_Ghosts = in->d_Ghosts; }
inline std::ostream &operator<<( std::ostream &out, const Vector::shared_ptr p ) { return operator<<( out, *p ); }
// clang-format on

/****************************************************************
 * Math API for Vector                                          *
 ****************************************************************/

inline void Vector::copy( const Vector &x ) { d_VectorOps->copy( *(x.getVectorData()), *(getVectorData()) ); }

inline void Vector::zero( void ) { d_VectorOps->zero( *(getVectorData()) ); }

inline void Vector::setToScalar( double alpha )
{
    d_VectorOps->setToScalar( alpha, *(getVectorData()) );
}

inline void Vector::setRandomValues( void ) { d_VectorOps->setRandomValues( *(getVectorData()) ); }

inline void Vector::setRandomValues( RNG::shared_ptr rng )
{
    d_VectorOps->setRandomValues( rng, *(getVectorData()) );
}

inline void Vector::scale( double alpha, const Vector &x )
{
    d_VectorOps->scale( alpha, *(x.getVectorData()), *(getVectorData()) );
}

inline void Vector::scale( double alpha ) { d_VectorOps->scale( alpha, *(getVectorData()) ); }

inline void Vector::add( const Vector &x, const Vector &y )
{
    d_VectorOps->add( *(x.getVectorData()), *(y.getVectorData()), *(getVectorData()) );
}

inline void Vector::subtract( const Vector &x, const Vector &y )
{
    d_VectorOps->subtract( *(x.getVectorData()), *(y.getVectorData()), *(getVectorData()) );
}

inline void Vector::multiply( const Vector &x, const Vector &y )
{
    d_VectorOps->multiply( *(x.getVectorData()), *(y.getVectorData()), *(getVectorData()) );
}

inline void Vector::divide( const Vector &x, const Vector &y )
{
    d_VectorOps->divide( *(x.getVectorData()), *(y.getVectorData()), *(getVectorData()) );
}

inline void Vector::reciprocal( const Vector &x )
{
    d_VectorOps->reciprocal( *(x.getVectorData()), *(getVectorData()) );
}

inline void Vector::linearSum( double alpha, const Vector &x, double beta, const Vector &y )
{
    d_VectorOps->linearSum( alpha, *(x.getVectorData()), beta, *(y.getVectorData()), *(getVectorData()) );
}

inline void Vector::axpy( double alpha, const Vector &x, const Vector &y )
{
    d_VectorOps->axpy( alpha, *(x.getVectorData()), *(y.getVectorData()), *(getVectorData()) );
}

inline void Vector::axpby( double alpha, double beta, const Vector &x )
{
    d_VectorOps->axpby( alpha, beta, *(x.getVectorData()), *(getVectorData()) );
}

inline void Vector::abs( const Vector &x ) { d_VectorOps->abs( *(x.getVectorData()), *(getVectorData()) ); }

inline void Vector::addScalar( const Vector &x, double alpha_in )
{
    d_VectorOps->addScalar( *(x.getVectorData()), alpha_in, *(getVectorData()) );
}

inline double Vector::min( void ) const { return d_VectorOps->min( *(getVectorData()) ); }

inline double Vector::max( void ) const { return d_VectorOps->max( *(getVectorData()) ); }

inline double Vector::L1Norm( void ) const { return d_VectorOps->L1Norm( *(getVectorData()) ); }

inline double Vector::L2Norm( void ) const { return d_VectorOps->L2Norm( *(getVectorData()) ); }

inline double Vector::maxNorm( void ) const { return d_VectorOps->maxNorm( *(getVectorData()) ); }

inline double Vector::minQuotient( const Vector &x ) const
{
    return d_VectorOps->minQuotient( *(x.getVectorData()), *(getVectorData()) );
}

inline double Vector::wrmsNorm( const Vector &x, const Vector &y ) const
{
    return d_VectorOps->wrmsNorm( *(x.getVectorData()), *(y.getVectorData()) );
}

inline double
Vector::wrmsNormMask( const Vector &x, const Vector &mask, const Vector &y ) const
{
    return d_VectorOps->wrmsNormMask( *(x.getVectorData()), mask, *(y.getVectorData()) );
}

inline double Vector::dot( const Vector &x ) const
{
  return d_VectorOps->dot( *(getVectorData()), *(x.getVectorData()) );
}

inline bool Vector::equals( const Vector &a, double tol ) const
{
    return d_VectorOps->equals( a, *(getVectorData()), tol );
}

inline double Vector::localMin( void ) const { return d_VectorOps->localMin( *(getVectorData()) ); }

inline double Vector::localMax( void ) const { return d_VectorOps->localMax( *(getVectorData()) ); }

inline double Vector::localL1Norm( void ) const
{
    return d_VectorOps->localL1Norm( *(getVectorData()) );
}

inline double Vector::localL2Norm( void ) const
{
    return d_VectorOps->localL2Norm( *(getVectorData()) );
}

inline double Vector::localMaxNorm( void ) const
{
    return d_VectorOps->localMaxNorm( *(getVectorData()) );
}

inline double Vector::localDot( const Vector &x ) const
{
    return d_VectorOps->localDot( *(x.getVectorData()), *(getVectorData()) );
}

inline double Vector::localMinQuotient( const Vector &x ) const
{
    return d_VectorOps->localMinQuotient( *(x.getVectorData()), *(getVectorData()) );
}

inline double Vector::localWrmsNorm( const Vector &x ) const
{
    return d_VectorOps->localWrmsNorm( *(x.getVectorData()), *(getVectorData()) );
}

inline double
Vector::localWrmsNormMask( const Vector &x, const Vector &mask, const Vector &y ) const
{
  return d_VectorOps->localWrmsNormMask( *(x.getVectorData()), *(mask.getVectorData()), *(y.getVectorData()) );
}

inline bool Vector::localEquals( const Vector &x, double tol ) const
{
    return d_VectorOps->localEquals( *(x.getVectorData()), *(getVectorData()), tol );
}

/****************************************************************
 * Wrappers for shared_ptr                                       *
 ****************************************************************/

inline bool Vector::equals( std::shared_ptr<const Vector> x, double tol ) const
{
    return equals( *x, tol );
}

inline void Vector::copy( std::shared_ptr<const Vector> x ) { copy( *x ); }

inline void Vector::scale( double alpha, std::shared_ptr<const Vector> x )
{
    scale( alpha, *x );
}

inline void Vector::add( std::shared_ptr<const Vector> x, std::shared_ptr<const Vector> y )
{
    add( *x, *y );
}

inline void Vector::addScalar( std::shared_ptr<const Vector> x, double alpha )
{
    addScalar( *x, alpha );
}

inline void Vector::subtract( std::shared_ptr<const Vector> x,
                              std::shared_ptr<const Vector> y )
{
    subtract( *x, *y );
}
inline void Vector::multiply( std::shared_ptr<const Vector> x,
                              std::shared_ptr<const Vector> y )
{
    multiply( *x, *y );
}
inline void Vector::divide( std::shared_ptr<const Vector> x,
                            std::shared_ptr<const Vector> y )
{
    divide( *x, *y );
}
inline void Vector::reciprocal( std::shared_ptr<const Vector> x ) { reciprocal( *x ); }
inline void Vector::linearSum( double alpha,
                               std::shared_ptr<const Vector> x,
                               double beta,
                               std::shared_ptr<const Vector> y )
{
    linearSum( alpha, *x, beta, *y );
}
inline void Vector::axpy( double alpha,
                          std::shared_ptr<const Vector> x,
                          std::shared_ptr<const Vector> y )
{
    axpy( alpha, *x, *y );
}
inline void Vector::axpby( double alpha, double beta, std::shared_ptr<const Vector> x )
{
    axpby( alpha, beta, *x );
}
inline void Vector::abs( std::shared_ptr<const Vector> x ) { return abs( *x ); }
inline double Vector::dot( std::shared_ptr<const Vector> x ) const { return dot( *x ); }

inline double Vector::minQuotient( std::shared_ptr<const Vector> x ) const
{
    return minQuotient( *x );
}
inline double Vector::wrmsNorm( std::shared_ptr<const Vector> x,
                                std::shared_ptr<const Vector> y ) const
{
    return d_VectorOps->wrmsNorm( *x, *y );
}
inline double Vector::wrmsNormMask( std::shared_ptr<const Vector> x,
                                    std::shared_ptr<const Vector> mask,
                                    std::shared_ptr<const Vector> y ) const
{
  return d_VectorOps->wrmsNormMask( *x, *mask, *y );
}

} // namespace LinearAlgebra
} // namespace AMP

#endif
