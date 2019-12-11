#include "AMP/utils/Utilities.h"
#include "VectorSelector.h"

namespace AMP {
namespace LinearAlgebra {


inline std::shared_ptr<VectorEngine> ManagedVector::getVectorEngine() { return d_Engine; }
inline std::shared_ptr<const VectorEngine> ManagedVector::getVectorEngine() const
{
    return d_Engine;
}


inline std::string ManagedVector::type() const
{
    if ( d_vBuffer )
        return " ( managed data )";
    std::string retVal = " ( managed view of ";
    auto vec           = std::dynamic_pointer_cast<Vector>( d_Engine );
    retVal += vec->type();
    retVal += " )";
    return retVal;
}


inline Vector::shared_ptr ManagedVector::getRootVector()
{
    if ( d_vBuffer )
        return shared_from_this();
    auto vec = std::dynamic_pointer_cast<ManagedVector>( d_Engine );
    if ( vec != nullptr )
        return vec->getRootVector();
    return std::dynamic_pointer_cast<Vector>( d_Engine )->shared_from_this();
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
        result = std::dynamic_pointer_cast<Vector>( d_Engine )->selectInto( s );
    }
    return result;
}


inline Vector::const_shared_ptr ManagedVector::selectInto( const VectorSelector &s ) const
{
    Vector::const_shared_ptr result;
    if ( d_vBuffer ) {
        result = Vector::selectInto( s );
    } else {
        result = std::dynamic_pointer_cast<Vector>( d_Engine )->selectInto( s );
    }
    return result;
}


inline ManagedVectorParameters::ManagedVectorParameters()
    : d_Buffer( nullptr ), d_CloneEngine( true )
{
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


inline std::shared_ptr<ParameterBase> ManagedVector::getParameters()
{
    return std::dynamic_pointer_cast<ParameterBase>( d_pParameters );
}


inline std::shared_ptr<ManagedVectorParameters> ManagedVector::getManagedVectorParameters()
{
    return d_pParameters;
}


inline size_t ManagedVector::getLocalSize() const { return d_Engine->getLocalSize(); }


inline size_t ManagedVector::getGlobalSize() const { return d_Engine->getGlobalSize(); }


inline ManagedVector::~ManagedVector() {}
} // namespace LinearAlgebra
} // namespace AMP
