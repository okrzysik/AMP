#ifndef included_AMP_Vector_inline
#define included_AMP_Vector_inline

#include "vectors/VectorDataIterator.h"
#include <algorithm>


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Get basic info                                                *
****************************************************************/
inline AMP::shared_ptr<ParameterBase> Vector::getParameters()
{
    return AMP::shared_ptr<ParameterBase>();
}
inline AMP::Discretization::DOFManager::shared_ptr Vector::getDOFManager() const
{
    return d_DOFManager;
}
inline AMP_MPI Vector::getComm() const { return d_CommList->getComm(); }


/****************************************************************
* Subset for variable name                                      *
****************************************************************/
inline Vector::shared_ptr Vector::subsetVectorForVariable( const std::string& name )
{
    auto var = AMP::make_shared<AMP::LinearAlgebra::Variable>(name);
    return subsetVectorForVariable( var );
}
inline Vector::const_shared_ptr
Vector::constSubsetVectorForVariable( const std::string& name ) const
{
    auto var = AMP::make_shared<AMP::LinearAlgebra::Variable>(name);
    return constSubsetVectorForVariable( var );
}


/****************************************************************
* getView/hasView                                               *
****************************************************************/
template <typename VIEW_TYPE>
Vector::shared_ptr Vector::getView() const
{
    for ( size_t i = 0; i != d_Views->size(); i++ ) {
        if ( ( *d_Views )[i].lock() ) {
            if ( ( *d_Views )[i].lock()->isA<VIEW_TYPE>() ) {
                return Vector::shared_ptr( ( *d_Views )[i] );
            }
        }
    }
    return Vector::shared_ptr();
}
template <typename VIEW_TYPE>
bool Vector::hasView() const
{
    for ( size_t i = 0; i != d_Views->size(); i++ ) {
        if ( ( *d_Views )[i].lock() ) {
            if ( ( *d_Views )[i].lock()->isA<VIEW_TYPE>() ) {
                return true;
            }
        }
    }
    return false;
}


/****************************************************************
* Misc functions                                                *
****************************************************************/
inline void Vector::setDefaultRNG( RNG::shared_ptr p ) { d_DefaultRNG = p; }

inline RNG::shared_ptr Vector::getDefaultRNG()
{
    if ( !d_DefaultRNG ) {
        AMP_MPI globalComm( AMP_COMM_WORLD );
        int rank = globalComm.getRank();
        RNGParameters::shared_ptr params(
            new RNGParameters( RNGParameters::RNGOptions::USE_GLOBAL_SEED, static_cast<size_t>( rank ) ) );
        d_DefaultRNG = RNG::shared_ptr( new RNG( params ) );
    }
    return d_DefaultRNG;
}

inline void Vector::requireSameSize( Vector &rhs )
{
    if ( rhs.getLocalSize() != getLocalSize() ) {
        AMP_ERROR( "Vectors are not of compatible size" );
    }
    if ( rhs.getGlobalSize() != getGlobalSize() ) {
        AMP_ERROR( "Vectors are not of compatible size" );
    }
}

inline const Variable::shared_ptr Vector::getVariable() const { return d_pVariable; }

inline Variable::shared_ptr Vector::getVariable()
{
    return d_pVariable; // Fix this!
}

inline Vector::shared_ptr Vector::cloneVector() const { return cloneVector( getVariable() ); }

inline void Vector::setVariable( const Variable::shared_ptr name )
{
    AMP_ASSERT( name.get() != nullptr );
    d_pVariable = name;
}

inline void Vector::addScalar( const_shared_ptr x, double alpha )
{
    Vector::shared_ptr one_vec = cloneVector();
    one_vec->setToScalar( 1. );
    axpy( alpha, one_vec, x );
}


/****************************************************************
* Wrappers for shared_ptr                                       *
****************************************************************/
// clang-format off
inline  bool Vector::equals( Vector::const_shared_ptr rhs, double tol ) const { return equals( *rhs, tol ); }
inline  void Vector::swapVectors( shared_ptr other ) { swapVectors( *other ); }
inline  void Vector::aliasVector( shared_ptr other ) { aliasVector( *other ); }
inline  void Vector::scale( double alpha, const_shared_ptr x ) { scale( alpha, *x ); }
inline  void Vector::add( const_shared_ptr x, const_shared_ptr y ) { add( *x, *y ); }
inline  void Vector::subtract( const_shared_ptr x, const_shared_ptr y ) { subtract( *x, *y ); }
inline  void Vector::multiply( const_shared_ptr x, const_shared_ptr y ) { multiply( *x, *y ); }
inline  void Vector::divide( const_shared_ptr x, const_shared_ptr y ) { divide( *x, *y ); }
inline  void Vector::reciprocal( const_shared_ptr x ) { reciprocal( *x ); }
inline double Vector::minQuotient( const_shared_ptr x, const_shared_ptr y ) { return minQuotient( *x, *y ); }
inline double Vector::wrmsNorm( const_shared_ptr x, const_shared_ptr y ) { return wrmsNorm( *x, *y ); }
inline  void Vector::linearSum( double alpha, const_shared_ptr x, double beta, const_shared_ptr y ) { linearSum( alpha, *x, beta, *y ); }
inline  void Vector::axpy( double alpha, const_shared_ptr x, const_shared_ptr y ) { axpy( alpha, *x, *y ); }
inline  void Vector::axpby( double alpha, double beta, const_shared_ptr x ) { axpby( alpha, beta, *x ); }
inline  void Vector::abs( const_shared_ptr x ) { this->abs( *x ); }
inline double Vector::dot( const_shared_ptr x ) const { return dot( *x ); }
inline  void Vector::addCommunicationListToParameters( CommunicationList::shared_ptr ) {}
inline  void Vector::aliasGhostBuffer( shared_ptr in ) { d_Ghosts = in->d_Ghosts; }
inline std::ostream &operator<<( std::ostream &out, const Vector::shared_ptr p ) { return operator<<( out, *p ); }
// clang-format on



} // LinearAlgebra namespace
} // AMP namespace

#endif
