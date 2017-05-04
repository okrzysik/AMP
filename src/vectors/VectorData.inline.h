#ifndef included_AMP_VectorData_inline
#define included_AMP_VectorData_inline

#include "vectors/VectorDataIterator.h"

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
* Create vector iterators                                       *
****************************************************************/
inline VectorDataIterator VectorData::begin() { return VectorDataIterator( this, 0 ); }
inline VectorDataIterator VectorData::end() { return VectorDataIterator( this, getLocalSize() ); }
inline ConstVectorDataIterator VectorData::begin() const { return ConstVectorDataIterator( this, 0 ); }
inline ConstVectorDataIterator VectorData::end() const { return ConstVectorDataIterator( this, getLocalSize() ); }
inline size_t VectorData::getGhostSize() const { return d_Ghosts->size(); }


/****************************************************************
* Update status                                                 *
****************************************************************/
inline VectorData::UpdateState VectorData::getUpdateStatus() const { return *d_UpdateState; }
inline void VectorData::setUpdateStatus( UpdateState state ) { *d_UpdateState = state; }
inline void VectorData::setUpdateStatusPtr( AMP::shared_ptr<UpdateState> rhs ) { d_UpdateState = rhs; }
inline AMP::shared_ptr<VectorData::UpdateState> VectorData::getUpdateStatusPtr() const
{
    return d_UpdateState;
}
inline void VectorData::dataChanged()
{
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}


/****************************************************************
* Templated functions                                           *
****************************************************************/
template <typename TYPE>
TYPE *VectorData::getRawDataBlock( size_t i )
{
    return static_cast<TYPE *>( this->getRawDataBlockAsVoid( i ) );
}
template <typename TYPE>
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


} // LinearAlgebra namespace
} // AMP namespace

#endif
