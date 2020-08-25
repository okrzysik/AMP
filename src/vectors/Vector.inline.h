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
inline void Vector::aliasGhostBuffer( shared_ptr in ) { d_VectorData->aliasGhostBuffer(in->d_VectorData); }
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
  return d_VectorOps->wrmsNormMask( *(x.getVectorData()), *(mask.getVectorData()), *(y.getVectorData()) );
}

inline double Vector::dot( const Vector &x ) const
{
  return d_VectorOps->dot( *(getVectorData()), *(x.getVectorData()) );
}

inline bool Vector::equals( const Vector &a, double tol ) const
{
  return d_VectorOps->equals( *(a.getVectorData()), *(getVectorData()), tol );
}

/****************************************************************
 * Get individual values                                     *
 ****************************************************************/
inline double Vector::getValueByGlobalID( size_t i ) const
{
    double ans;
    getValuesByGlobalID( 1, &i, &ans );
    return ans;
}
inline double Vector::getLocalValueByGlobalID( size_t i ) const
{
    double ans;
    getLocalValuesByGlobalID( 1, &i, &ans );
    return ans;
}
inline double Vector::getGhostValueByGlobalID( size_t i ) const
{
    double ans;
    getGhostValuesByGlobalID( 1, &i, &ans );
    return ans;
}
inline double Vector::getValueByLocalID( size_t ndx ) const
{
    double ans;
    getValuesByLocalID( 1, &ndx, &ans );
    return ans;
}

} // namespace LinearAlgebra
} // namespace AMP

#endif
