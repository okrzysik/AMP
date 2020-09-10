#include "AMP/vectors/ManagedVector.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <typeinfo>

#include "AMP/utils/Utilities.h"
#include "AMP/vectors/ManagedVectorData.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/operations/ManagedVectorOperations.h"

namespace AMP {
namespace LinearAlgebra {


// Helper functions
static inline ManagedVector *getManaged( Vector *x )
{
    auto y = dynamic_cast<ManagedVector *>( x );
    AMP_INSIST( y != nullptr, "x is not a ManagedVector" );
    return y;
}
static inline std::shared_ptr<ManagedVector> getManaged( std::shared_ptr<Vector> x )
{
    auto y = std::dynamic_pointer_cast<ManagedVector>( x );
    AMP_INSIST( y != nullptr, "x is not a ManagedVector" );
    return y;
}


/********************************************************
 * Constructors                                          *
 ********************************************************/
ManagedVector::ManagedVector( VectorParameters::shared_ptr params_in )
    : Vector( params_in ),
      d_pParameters( std::dynamic_pointer_cast<ManagedVectorParameters>( params_in ) )
{
    d_VectorOps = std::make_shared<ManagedVectorOperations>();
    setVectorData( std::make_shared<ManagedVectorData>( params_in ) );
}
ManagedVector::ManagedVector( shared_ptr alias )
    : Vector( std::dynamic_pointer_cast<VectorParameters>( getManaged( alias )->getParameters() ) )
{
    auto vec = getManaged( alias );
    setVectorData( vec->d_VectorData );
    d_VectorOps   = vec->d_VectorOps;
    d_pParameters = vec->d_pParameters;
    setVariable( vec->getVariable() );
}
ManagedVector::~ManagedVector() {}


/********************************************************
 * Subset                                                *
 ********************************************************/
Vector::shared_ptr ManagedVector::subsetVectorForVariable( Variable::const_shared_ptr name )
{
    Vector::shared_ptr retVal;
    if ( !retVal )
        retVal = Vector::subsetVectorForVariable( name );
    if ( !retVal ) {
        auto vec = getVectorEngine();
        if ( vec )
            retVal = vec->subsetVectorForVariable( name );
    }
    return retVal;
}
Vector::const_shared_ptr
ManagedVector::constSubsetVectorForVariable( Variable::const_shared_ptr name ) const
{
    Vector::const_shared_ptr retVal;
    if ( !retVal )
        retVal = Vector::constSubsetVectorForVariable( name );
    if ( !retVal ) {
        auto const vec = getVectorEngine();
        if ( vec )
            retVal = vec->constSubsetVectorForVariable( name );
    }
    return retVal;
}


bool ManagedVector::isAnAliasOf( Vector &rhs )
{
    bool retVal = false;

    auto other = getManaged( &rhs );
    if ( other != nullptr ) {
        retVal = std::dynamic_pointer_cast<ManagedVectorData>( d_VectorData )
                     ->isAnAliasOf( *( rhs.getVectorData() ) );
    }
    return retVal;
}

bool ManagedVector::isAnAliasOf( Vector::shared_ptr rhs ) { return isAnAliasOf( *rhs ); }

void ManagedVector::swapVectors( Vector &other )
{
    d_VectorData->swapData( *other.getVectorData() );
    auto in = getManaged( &other );
    std::swap( d_pParameters, in->d_pParameters );
}

std::shared_ptr<Vector> ManagedVector::cloneVector( const Variable::shared_ptr name ) const
{
    std::shared_ptr<ManagedVector> retVal( getNewRawPtr() );
    retVal->setVectorData( d_VectorData->cloneData() );
    retVal->d_DOFManager = getDOFManager();
    retVal->setVariable( name );
    return retVal;
}

std::string ManagedVector::type() const
{
    AMP_ASSERT( d_VectorData );
    return d_VectorData->VectorDataName();
}

Vector::shared_ptr ManagedVector::getRootVector()
{
    if ( std::dynamic_pointer_cast<ManagedVectorData>( d_VectorData )->hasBuffer() )
        return shared_from_this();

    auto engine = getVectorEngine();
    auto vec    = std::dynamic_pointer_cast<ManagedVector>( engine );
    if ( vec != nullptr ) {
        auto rvec = vec->getRootVector();
        AMP_INSIST( rvec->getCommunicationList(),
                    "Root vector does not have a communication list" );
    }

    AMP_INSIST( engine->getCommunicationList(),
                "Managed vector engine does not have a communication list" );
    return engine->shared_from_this();
}

Vector::shared_ptr ManagedVector::selectInto( const VectorSelector &s )
{
    return Vector::selectInto( s );
}

Vector::const_shared_ptr ManagedVector::selectInto( const VectorSelector &s ) const
{
    return Vector::selectInto( s );
}

std::shared_ptr<ParameterBase> ManagedVector::getParameters()
{
    return std::dynamic_pointer_cast<ParameterBase>( d_pParameters );
}

Vector::shared_ptr ManagedVector::getVectorEngine( void )
{
    AMP_ASSERT( d_VectorData );
    auto data = std::dynamic_pointer_cast<ManagedVectorData>( d_VectorData );
    AMP_ASSERT( data );
    return data->getVectorEngine();
}

Vector::const_shared_ptr ManagedVector::getVectorEngine( void ) const
{
    AMP_ASSERT( d_VectorData );
    const auto data = std::dynamic_pointer_cast<const ManagedVectorData>( d_VectorData );
    AMP_ASSERT( data );
    return data->getVectorEngine();
}

} // namespace LinearAlgebra
} // namespace AMP
