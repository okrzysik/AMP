#include "VectorSelector.h"
#include "utils/Utilities.h"

namespace AMP {
namespace LinearAlgebra {


inline VectorEngine::shared_ptr ManagedVector::getVectorEngine() { return d_Engine; }
inline VectorEngine::const_shared_ptr ManagedVector::getVectorEngine() const { return d_Engine; }


inline std::string ManagedVector::type() const
{
    if ( d_vBuffer )
        return " ( managed data )";
    std::string retVal = " ( managed view of ";
    retVal += d_Engine->castTo<Vector>().type();
    retVal += " )";
    return retVal;
}


inline Vector::shared_ptr ManagedVector::getRootVector()
{
    if ( d_vBuffer )
        return shared_from_this();
    Vector *t = &( d_Engine->castTo<Vector>() );
    if ( t->isA<ManagedVector>() )
        return t->castTo<ManagedVector>().getRootVector();
    return t->shared_from_this();
}


inline Vector::const_iterator ManagedVector::begin() const
{
    if ( d_vBuffer )
        return Vector::begin();
    else
        return d_Engine->castTo<const Vector>().begin();
}


inline Vector::const_iterator ManagedVector::end() const
{
    if ( d_vBuffer )
        return Vector::end();
    else
        return d_Engine->castTo<const Vector>().end();
}


inline Vector::iterator ManagedVector::begin()
{
    if ( d_vBuffer )
        return Vector::begin();
    else
        return d_Engine->castTo<Vector>().begin();
}


inline Vector::iterator ManagedVector::end()
{
    if ( d_vBuffer )
        return Vector::end();
    else
        return d_Engine->castTo<Vector>().end();
}


inline void ManagedVector::dataChanged()
{
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}


inline Vector::shared_ptr ManagedVector::selectInto( const VectorSelector &s )
{
    Vector::shared_ptr result;
    if ( d_vBuffer ) {
        result = Vector::selectInto( s );
    } else {
        result = d_Engine->castTo<Vector>().selectInto( s );
    }
    return result;
}


inline Vector::const_shared_ptr ManagedVector::selectInto( const VectorSelector &s ) const
{
    Vector::const_shared_ptr result;
    if ( d_vBuffer ) {
        result = Vector::selectInto( s );
    } else {
        result = d_Engine->castTo<Vector>().selectInto( s );
    }
    return result;
}


inline ManagedVectorParameters::ManagedVectorParameters()
{
    d_CloneEngine = true;
    d_Buffer      = VectorEngine::BufferPtr();
}


inline void *ManagedVector::getRawDataBlockAsVoid( size_t i )
{
    return d_Engine->getDataBlock( i );
}


inline const void *ManagedVector::getRawDataBlockAsVoid( size_t i ) const
{
    return d_Engine->getDataBlock( i );
}


inline void ManagedVector::addCommunicationListToParameters( CommunicationList::shared_ptr comm )
{
    d_pParameters->d_CommList = comm;
}


inline size_t ManagedVector::numberOfDataBlocks() const { return d_Engine->numberOfDataBlocks(); }


inline size_t ManagedVector::sizeOfDataBlock( size_t i ) const
{
    return d_Engine->sizeOfDataBlock( i );
}


inline bool ManagedVector::isAnAliasOf( Vector::shared_ptr rhs ) { return isAnAliasOf( *rhs ); }


inline AMP::shared_ptr<ParameterBase> ManagedVector::getParameters()
{
    return AMP::dynamic_pointer_cast<ParameterBase>( d_pParameters );
}


inline AMP::shared_ptr<ManagedVectorParameters> ManagedVector::getManagedVectorParameters()
{
    return d_pParameters;
}


inline size_t ManagedVector::getLocalSize() const { return d_Engine->getLocalSize(); }


inline size_t ManagedVector::getGlobalSize() const { return d_Engine->getGlobalSize(); }


inline AMP::shared_ptr<Vector> ManagedVector::cloneVector( const Variable::shared_ptr name ) const
{
    AMP::shared_ptr<Vector> retVal( getNewRawPtr() );
    if ( !d_vBuffer ) {
        retVal->castTo<ManagedVector>().d_Engine =
            d_Engine->cloneEngine( VectorEngine::BufferPtr() );
    }
    retVal->setVariable( name );
    return retVal;
}


inline ManagedVector::~ManagedVector() {}
}
}
