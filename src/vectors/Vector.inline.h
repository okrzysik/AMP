#ifndef included_AMP_Vector_inline
#define included_AMP_Vector_inline

#include "AMP/vectors/data/VectorDataIterator.h"
#include <algorithm>


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Get basic info                                                *
 ****************************************************************/
inline std::shared_ptr<AMP::Discretization::DOFManager> Vector::getDOFManager() const
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
std::shared_ptr<VIEW_TYPE> Vector::getView() const
{
    typedef typename std::remove_cv<VIEW_TYPE>::type TYPE;
    for ( size_t i = 0; i != d_Views->size(); i++ ) {
        if ( ( *d_Views )[i].type() == typeid( std::weak_ptr<TYPE> ) ) {
            auto ptr = std::any_cast<std::weak_ptr<TYPE>>( ( *d_Views )[i] );
            auto vec = ptr.lock();
            if ( vec )
                return vec;
        }
    }
    return std::shared_ptr<VIEW_TYPE>();
}
template<typename VIEW_TYPE>
bool Vector::hasView() const
{
    return getView<VIEW_TYPE>() != nullptr;
}
template<typename VIEW_TYPE>
void Vector::registerView( std::shared_ptr<VIEW_TYPE> v ) const
{
    typedef typename std::remove_cv<VIEW_TYPE>::type TYPE;
    for ( size_t i = 0; i != d_Views->size(); i++ ) {
        if ( ( *d_Views )[i].type() == typeid( std::weak_ptr<TYPE> ) ) {
            auto ptr = std::any_cast<std::weak_ptr<TYPE>>( ( *d_Views )[i] );
            auto vec = ptr.lock();
            if ( vec == v )
                return;
        }
    }
    std::weak_ptr<TYPE> ptr = v;
    d_Views->push_back( std::any( ptr ) );
}


/****************************************************************
 * Misc functions                                                *
 ****************************************************************/
inline const std::shared_ptr<Variable> Vector::getVariable() const { return d_pVariable; }
inline std::shared_ptr<Variable> Vector::getVariable() { return d_pVariable; }
inline void Vector::setVariable( const std::shared_ptr<Variable> name )
{
    AMP_ASSERT( name );
    d_pVariable = name;
}


/****************************************************************
 * Wrappers for shared_ptr                                       *
 ****************************************************************/
// clang-format off
inline void Vector::swapVectors( shared_ptr other ) { swapVectors( *other ); }
inline void Vector::aliasGhostBuffer( shared_ptr in ) { d_VectorData->aliasGhostBuffer(in->d_VectorData); }
inline std::ostream &operator<<( std::ostream &out, const Vector::shared_ptr p ) { return operator<<( out, *p ); }
// clang-format on


/****************************************************************
 * Math API for Vector                                          *
 ****************************************************************/
inline void Vector::copy( const Vector &x )
{
    d_VectorOps->copy( *x.getVectorData(), *getVectorData() );
}
inline void Vector::zero( void ) { d_VectorOps->zero( *getVectorData() ); }
inline void Vector::setToScalar( const Scalar &alpha )
{
    d_VectorOps->setToScalar( alpha, *getVectorData() );
}
inline void Vector::scale( const Scalar &alpha, const Vector &x )
{
    d_VectorOps->scale( alpha, *x.getVectorData(), *getVectorData() );
}
inline void Vector::scale( const Scalar &alpha ) { d_VectorOps->scale( alpha, *getVectorData() ); }
inline void Vector::add( const Vector &x, const Vector &y )
{
    d_VectorOps->add( *x.getVectorData(), *y.getVectorData(), *getVectorData() );
}
inline void Vector::subtract( const Vector &x, const Vector &y )
{
    d_VectorOps->subtract( *x.getVectorData(), *y.getVectorData(), *getVectorData() );
}
inline void Vector::multiply( const Vector &x, const Vector &y )
{
    d_VectorOps->multiply( *x.getVectorData(), *y.getVectorData(), *getVectorData() );
}
inline void Vector::divide( const Vector &x, const Vector &y )
{
    d_VectorOps->divide( *x.getVectorData(), *y.getVectorData(), *getVectorData() );
}
inline void Vector::reciprocal( const Vector &x )
{
    d_VectorOps->reciprocal( *x.getVectorData(), *getVectorData() );
}
inline void
Vector::linearSum( const Scalar &alpha, const Vector &x, const Scalar &beta, const Vector &y )
{
    d_VectorOps->linearSum( alpha, *x.getVectorData(), beta, *y.getVectorData(), *getVectorData() );
}
inline void Vector::axpy( const Scalar &alpha, const Vector &x, const Vector &y )
{
    d_VectorOps->axpy( alpha, *x.getVectorData(), *y.getVectorData(), *getVectorData() );
}
inline void Vector::axpby( const Scalar &alpha, const Scalar &beta, const Vector &x )
{
    d_VectorOps->axpby( alpha, beta, *x.getVectorData(), *getVectorData() );
}
inline void Vector::abs( const Vector &x )
{
    d_VectorOps->abs( *x.getVectorData(), *getVectorData() );
}
inline void Vector::addScalar( const Vector &x, const Scalar &alpha_in )
{
    d_VectorOps->addScalar( *x.getVectorData(), alpha_in, *getVectorData() );
}
inline Scalar Vector::min( void ) const { return d_VectorOps->min( *getVectorData() ); }
inline Scalar Vector::max( void ) const { return d_VectorOps->max( *getVectorData() ); }
inline Scalar Vector::L1Norm( void ) const { return d_VectorOps->L1Norm( *getVectorData() ); }
inline Scalar Vector::L2Norm( void ) const { return d_VectorOps->L2Norm( *getVectorData() ); }
inline Scalar Vector::maxNorm( void ) const { return d_VectorOps->maxNorm( *getVectorData() ); }
inline Scalar Vector::minQuotient( const Vector &x ) const
{
    return d_VectorOps->minQuotient( *x.getVectorData(), *getVectorData() );
}
inline Scalar Vector::wrmsNorm( const Vector &x, const Vector &y ) const
{
    return d_VectorOps->wrmsNorm( *x.getVectorData(), *y.getVectorData() );
}
inline Scalar Vector::wrmsNormMask( const Vector &x, const Vector &mask, const Vector &y ) const
{
    return d_VectorOps->wrmsNormMask(
        *x.getVectorData(), *mask.getVectorData(), *y.getVectorData() );
}
inline Scalar Vector::dot( const Vector &x ) const
{
    return d_VectorOps->dot( *getVectorData(), *x.getVectorData() );
}
inline bool Vector::equals( const Vector &a, const Scalar &tol ) const
{
    return d_VectorOps->equals( *a.getVectorData(), *getVectorData(), tol );
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
