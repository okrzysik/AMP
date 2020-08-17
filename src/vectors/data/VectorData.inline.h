#ifndef included_AMP_VectorData_inline
#define included_AMP_VectorData_inline

#include "AMP/vectors/data/VectorDataIterator.h"

#include <algorithm>


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Get the size of the vector                                    *
 ****************************************************************/
inline size_t VectorData::getGlobalMaxID() const { return getGlobalSize(); }
inline size_t VectorData::getLocalMaxID() const { return getLocalSize(); }
inline size_t VectorData::getLocalStartID() const { return d_CommList->getStartGID(); }
inline CommunicationList::shared_ptr VectorData::getCommunicationList() const { return d_CommList; }
inline void VectorData::aliasGhostBuffer( std::shared_ptr<VectorData> in ) { d_Ghosts = in->d_Ghosts; }
inline bool VectorData::containsGlobalElement( size_t i )
{
    if ( ( i >= d_CommList->getStartGID() ) &&
         ( i < d_CommList->getStartGID() + d_CommList->numLocalRows() ) )
        return true;
    return std::find( d_CommList->getGhostIDList().begin(),
                      d_CommList->getGhostIDList().end(),
                      i ) != d_CommList->getGhostIDList().end();
}


/****************************************************************
 * Get the type of data                                          *
 ****************************************************************/
template<typename TYPE>
bool VectorData::isType() const
{
    bool test = true;
    auto hash = typeid( TYPE ).hash_code();
    for ( size_t i = 0; i < numberOfDataBlocks(); i++ )
        test = test && isTypeId( hash, i );
    return test;
}
template<typename TYPE>
bool VectorData::isBlockType( size_t i ) const
{
    auto hash = typeid( TYPE ).hash_code();
    return isTypeId( hash, i );
}


/****************************************************************
 * Create vector iterators                                       *
 ****************************************************************/
template<class TYPE>
inline VectorDataIterator<TYPE> VectorData::begin()
{
    dataChanged();
    return VectorDataIterator<TYPE>( this, 0 );
}
template<class TYPE>
inline VectorDataIterator<const TYPE> VectorData::begin() const
{
    return VectorDataIterator<const TYPE>( const_cast<VectorData *>( this ), 0 );
}
template<class TYPE>
inline VectorDataIterator<const TYPE> VectorData::constBegin() const
{
    return VectorDataIterator<const TYPE>( const_cast<VectorData *>( this ), 0 );
}
template<class TYPE>
inline VectorDataIterator<TYPE> VectorData::end()
{
    dataChanged();
    return VectorDataIterator<TYPE>( this, getLocalSize() );
}
template<class TYPE>
inline VectorDataIterator<const TYPE> VectorData::constEnd() const
{
    return VectorDataIterator<const TYPE>( const_cast<VectorData *>( this ), getLocalSize() );
}
template<class TYPE>
inline VectorDataIterator<const TYPE> VectorData::end() const
{
    return VectorDataIterator<const TYPE>( const_cast<VectorData *>( this ), getLocalSize() );
}
inline size_t VectorData::getGhostSize() const { return d_Ghosts->size(); }


/****************************************************************
 * Update status                                                 *
 ****************************************************************/
inline VectorData::UpdateState VectorData::getUpdateStatus() const { return *d_UpdateState; }
inline void VectorData::setUpdateStatus( UpdateState state ) { *d_UpdateState = state; }
inline void VectorData::setUpdateStatusPtr( std::shared_ptr<UpdateState> rhs )
{
    d_UpdateState = rhs;
}
inline std::shared_ptr<VectorData::UpdateState> VectorData::getUpdateStatusPtr() const
{
    return d_UpdateState;
}


/****************************************************************
 * Templated functions                                           *
 ****************************************************************/
template<typename TYPE>
TYPE *VectorData::getRawDataBlock( size_t i )
{
    return static_cast<TYPE *>( this->getRawDataBlockAsVoid( i ) );
}
template<typename TYPE>
const TYPE *VectorData::getRawDataBlock( size_t i ) const
{
    return static_cast<const TYPE *>( this->getRawDataBlockAsVoid( i ) );
}


/****************************************************************
 * Set/Get individual values                                     *
 ****************************************************************/
inline void VectorData::setValueByLocalID( size_t i, const double val )
{
    setValuesByLocalID( 1, &i, &val );
}
inline void VectorData::setLocalValueByGlobalID( size_t i, const double val )
{
    setLocalValuesByGlobalID( 1, &i, &val );
}
inline void VectorData::setGhostValueByGlobalID( size_t i, const double val )
{
    setGhostValuesByGlobalID( 1, &i, &val );
}
inline void VectorData::setValueByGlobalID( size_t i, const double val )
{
    setValuesByGlobalID( 1, &i, &val );
}
inline void VectorData::addValueByLocalID( size_t i, const double val )
{
    addValuesByLocalID( 1, &i, &val );
}
inline void VectorData::addLocalValueByGlobalID( size_t i, const double val )
{
    addLocalValuesByGlobalID( 1, &i, &val );
}
inline void VectorData::addValueByGlobalID( size_t i, const double val )
{
    addValuesByGlobalID( 1, &i, &val );
}
inline double VectorData::getValueByGlobalID( size_t i ) const
{
    double ans;
    getValuesByGlobalID( 1, &i, &ans );
    return ans;
}
inline double VectorData::getLocalValueByGlobalID( size_t i ) const
{
    double ans;
    getLocalValuesByGlobalID( 1, &i, &ans );
    return ans;
}
inline double VectorData::getGhostValueByGlobalID( size_t i ) const
{
    double ans;
    getGhostValuesByGlobalID( 1, &i, &ans );
    return ans;
}
inline double VectorData::getValueByLocalID( size_t ndx ) const
{
    double ans;
    getValuesByLocalID( 1, &ndx, &ans );
    return ans;
}


} // namespace LinearAlgebra
} // namespace AMP

#endif
