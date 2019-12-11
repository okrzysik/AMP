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


} // namespace LinearAlgebra
} // namespace AMP

#endif
