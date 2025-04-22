#ifndef included_AMP_GhostDataHelper_hpp
#define included_AMP_GhostDataHelper_hpp

#include "AMP/IO/RestartManager.h"
#include "AMP/vectors/data/GhostDataHelper.h"

namespace AMP::LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
template<class TYPE, class Allocator>
GhostDataHelper<TYPE, Allocator>::GhostDataHelper()
    : d_UpdateState{ std::make_shared<UpdateState>() }
{
    *d_UpdateState = UpdateState::UNCHANGED;
}
template<class TYPE, class Allocator>
GhostDataHelper<TYPE, Allocator>::GhostDataHelper( std::shared_ptr<CommunicationList> list )
    : d_UpdateState{ std::make_shared<UpdateState>() }
{
    *d_UpdateState = UpdateState::UNCHANGED;
    setCommunicationList( list );
}


/****************************************************************
 * Get/Set the communication list                                *
 ****************************************************************/
template<class TYPE, class Allocator>
std::shared_ptr<CommunicationList> GhostDataHelper<TYPE, Allocator>::getCommunicationList() const
{
    return d_CommList;
}
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::setCommunicationList(
    std::shared_ptr<CommunicationList> comm )
{
    AMP_ASSERT( comm );
    d_CommList = comm;
    if ( d_CommList ) {
        d_Ghosts =
            std::make_shared<std::vector<double>>( d_CommList->getVectorReceiveBufferSize() );
        d_AddBuffer =
            std::make_shared<std::vector<double>>( d_CommList->getVectorReceiveBufferSize() );
    }
}


/****************************************************************
 * makeConsistent / UpdateState                                  *
 ****************************************************************/
template<class TYPE, class Allocator>
UpdateState GhostDataHelper<TYPE, Allocator>::getLocalUpdateStatus() const
{
    return *d_UpdateState;
}
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::setUpdateStatus( UpdateState state )
{
    *d_UpdateState = state;
}
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::setUpdateStatusPtr( std::shared_ptr<UpdateState> rhs )
{
    d_UpdateState = rhs;
}
template<class TYPE, class Allocator>
std::shared_ptr<UpdateState> GhostDataHelper<TYPE, Allocator>::getUpdateStatusPtr() const
{
    return d_UpdateState;
}
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::dataChanged()
{
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
    fireDataChange();
}
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::makeConsistent( ScatterType t )
{
    PROFILE( "makeConsistent" );
    if ( d_CommList ) {
        if ( t == ScatterType::CONSISTENT_ADD ) {
            AMP_ASSERT( *d_UpdateState != UpdateState::SETTING );
            d_CommList->scatter_add( *this );
            for ( auto &elem : *d_AddBuffer )
                elem = 0.0;
        }
        *d_UpdateState = UpdateState::SETTING;
        d_CommList->scatter_set( *this );
        *d_UpdateState = UpdateState::UNCHANGED;
    }
    this->setUpdateStatus( UpdateState::UNCHANGED );
}


/****************************************************************
 * Get ghost size                                                *
 ****************************************************************/
template<class TYPE, class Allocator>
size_t GhostDataHelper<TYPE, Allocator>::getGhostSize() const
{
    if ( d_Ghosts )
        return d_Ghosts->size();
    return 0;
}


/****************************************************************
 * Alias ghost buffer                                            *
 ****************************************************************/
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::aliasGhostBuffer( std::shared_ptr<VectorData> in )
{
    auto ghostData = std::dynamic_pointer_cast<GhostDataHelper<TYPE, Allocator>>( in );
    AMP_ASSERT( ghostData );
    d_Ghosts = ghostData->d_Ghosts;
}


/****************************************************************
 * Zero ghost buffers                                            *
 ****************************************************************/
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::zeroGhosts()
{
    if ( d_Ghosts ) {
        for ( auto &x : *d_Ghosts )
            x = 0;
    }
    if ( d_AddBuffer ) {
        for ( auto &x : *d_AddBuffer )
            x = 0;
    }
}


/****************************************************************
 * Check if vector contains a particular element                 *
 ****************************************************************/
template<class TYPE, class Allocator>
bool GhostDataHelper<TYPE, Allocator>::containsGlobalElement( size_t i ) const
{
    if ( ( i >= d_CommList->getStartGID() ) &&
         ( i < d_CommList->getStartGID() + d_CommList->numLocalRows() ) )
        return true;
    return std::find( d_CommList->getGhostIDList().begin(),
                      d_CommList->getGhostIDList().end(),
                      i ) != d_CommList->getGhostIDList().end();
}


/****************************************************************
 * Get/Set ghost values by global id                             *
 ****************************************************************/
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::setGhostValuesByGlobalID( size_t N,
                                                                 const size_t *ndx,
                                                                 const void *vals,
                                                                 const typeID &id )
{
    if ( id == AMP::getTypeID<double>() ) {
        auto data = reinterpret_cast<const double *>( vals );
        AMP_ASSERT( *d_UpdateState != UpdateState::ADDING );
        *d_UpdateState = UpdateState::SETTING;
        for ( size_t i = 0; i < N; i++ ) {
            if ( ( ndx[i] < d_localStart ) || ( ndx[i] >= ( d_localStart + d_localSize ) ) ) {
                ( *d_Ghosts )[d_CommList->getLocalGhostID( ndx[i] )] = data[i];
            } else {
                AMP_ERROR( "Non ghost index" );
            }
        }
    } else {
        AMP_ERROR( "Ghosts other than double are not supported yet" );
    }
}
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::addGhostValuesByGlobalID( size_t N,
                                                                 const size_t *ndx,
                                                                 const void *vals,
                                                                 const typeID &id )
{
    if ( id == AMP::getTypeID<double>() ) {
        auto data = reinterpret_cast<const double *>( vals );
        AMP_ASSERT( *d_UpdateState != UpdateState::SETTING );
        *d_UpdateState = UpdateState::ADDING;
        for ( size_t i = 0; i < N; i++ ) {
            if ( ( ndx[i] < d_localStart ) || ( ndx[i] >= ( d_localStart + d_localSize ) ) ) {
                ( *d_AddBuffer )[d_CommList->getLocalGhostID( ndx[i] )] += data[i];
            } else {
                AMP_ERROR( "Non ghost index" );
            }
        }
    } else {
        AMP_ERROR( "Ghosts other than double are not supported yet" );
    }
}
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::getGhostValuesByGlobalID( size_t N,
                                                                 const size_t *ndx,
                                                                 void *vals,
                                                                 const typeID &id ) const
{
    if ( id == AMP::getTypeID<double>() ) {
        auto data = reinterpret_cast<double *>( vals );
        for ( size_t i = 0; i < N; i++ ) {
            if ( ( ndx[i] < d_localStart ) || ( ndx[i] >= ( d_localStart + d_localSize ) ) ) {
                data[i] = ( *d_Ghosts )[d_CommList->getLocalGhostID( ndx[i] )] +
                          ( *d_AddBuffer )[d_CommList->getLocalGhostID( ndx[i] )];
            } else {
                AMP_ERROR( "Tried to get a non-ghost ghost value" );
            }
        }
    } else {
        AMP_ERROR( "Ghosts other than double are not supported yet" );
    }
}
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::getGhostAddValuesByGlobalID( size_t N,
                                                                    const size_t *ndx,
                                                                    void *vals,
                                                                    const typeID &id ) const
{
    if ( id == AMP::getTypeID<double>() ) {
        auto data = reinterpret_cast<double *>( vals );
        for ( size_t i = 0; i < N; i++ ) {
            if ( ( ndx[i] < d_localStart ) || ( ndx[i] >= ( d_localStart + d_localSize ) ) ) {
                data[i] = ( *d_AddBuffer )[d_CommList->getLocalGhostID( ndx[i] )];
            } else {
                AMP_ERROR( "Tried to get a non-ghost ghost value" );
            }
        }
    } else {
        AMP_ERROR( "Ghosts other than double are not supported yet" );
    }
}


/****************************************************************
 * dump data to ostream                                          *
 ****************************************************************/
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::dumpGhostedData( std::ostream &out, size_t offset ) const
{
    if ( !getCommunicationList() )
        return;
    const std::vector<size_t> &ghosts = getCommunicationList()->getGhostIDList();
    auto curVal                       = d_Ghosts->begin();
    for ( auto &ghost : ghosts ) {
        out << "  GID: " << ( ghost + offset ) << "  Value: " << ( *curVal ) << "\n";
        ++curVal;
    }
}


/****************************************************************
 * Default clone                                                 *
 ****************************************************************/
template<class TYPE, class Allocator>
const AMP_MPI &GhostDataHelper<TYPE, Allocator>::getComm() const
{
    AMP_ASSERT( d_CommList );
    return d_CommList->getComm();
}


/****************************************************************
 * dump data to ostream                                          *
 ****************************************************************/
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::copyGhostValues( const VectorData &rhs )
{
    if ( getGhostSize() == 0 ) {
        // No ghosts to fill, copy the consistency state from the rhs
        *d_UpdateState = rhs.getLocalUpdateStatus();
    } else if ( getGhostSize() == rhs.getGhostSize() ) {
        // The ghosts in the src vector match the current vector
        // Copy the ghosts from the rhs
        auto ghostIDs = getCommunicationList()->getGhostIDList();
        std::vector<double> values( ghostIDs.size() );
        rhs.getGhostValuesByGlobalID( ghostIDs.size(), &ghostIDs[0], &values[0] );
        this->setGhostValuesByGlobalID( ghostIDs.size(), &ghostIDs[0], &values[0] );
        // Copy the consistency state from the rhs
        *d_UpdateState = rhs.getLocalUpdateStatus();
    } else {
        // We can't copy the ghosts from the rhs
        // Use makeConsistent to fill the ghosts
        // Note: this will insure global communication
        *d_UpdateState = rhs.getLocalUpdateStatus();
        if ( *d_UpdateState == UpdateState::UNCHANGED )
            *d_UpdateState = UpdateState::LOCAL_CHANGED;
    }
}


/****************************************************************
 * Write/Read restart data                                       *
 ****************************************************************/
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::registerChildObjects(
    AMP::IO::RestartManager *manager ) const
{
    VectorData::registerChildObjects( manager );
    if ( d_CommList ) {
        auto id = manager->registerObject( d_CommList );
        AMP_ASSERT( id == d_CommList->getID() );
    }
    if ( d_UpdateState ) {
        auto id = manager->registerObject( d_UpdateState );
        AMP_ASSERT( id == reinterpret_cast<uint64_t>( d_UpdateState.get() ) );
    }
    if ( d_Ghosts ) {
        auto id = manager->registerObject( d_Ghosts );
        AMP_ASSERT( id == reinterpret_cast<uint64_t>( d_Ghosts.get() ) );
    }
    if ( d_AddBuffer ) {
        auto id = manager->registerObject( d_AddBuffer );
        AMP_ASSERT( id == reinterpret_cast<uint64_t>( d_AddBuffer.get() ) );
    }
}
template<class TYPE, class Allocator>
void GhostDataHelper<TYPE, Allocator>::writeRestart( int64_t fid ) const
{
    VectorData::writeRestart( fid );
    uint64_t commListID  = d_CommList ? d_CommList->getID() : 0;
    uint64_t updateID    = reinterpret_cast<uint64_t>( d_UpdateState.get() );
    uint64_t ghostID     = reinterpret_cast<uint64_t>( d_Ghosts.get() );
    uint64_t addBufferID = reinterpret_cast<uint64_t>( d_AddBuffer.get() );
    IO::writeHDF5( fid, "commListID", commListID );
    IO::writeHDF5( fid, "updateID", updateID );
    IO::writeHDF5( fid, "ghostID", ghostID );
    IO::writeHDF5( fid, "addBufferID", addBufferID );
}
template<class TYPE, class Allocator>
GhostDataHelper<TYPE, Allocator>::GhostDataHelper( int64_t fid, AMP::IO::RestartManager *manager )
    : VectorData( fid, manager )
{
    uint64_t commListID, updateID, ghostID, addBufferID;
    IO::readHDF5( fid, "commListID", commListID );
    IO::readHDF5( fid, "updateID", updateID );
    IO::readHDF5( fid, "ghostID", ghostID );
    IO::readHDF5( fid, "addBufferID", addBufferID );
    if ( commListID != 0 )
        d_CommList = manager->getData<CommunicationList>( commListID );
    if ( updateID != 0 )
        d_UpdateState = manager->getData<UpdateState>( updateID );
    if ( ghostID != 0 )
        d_Ghosts = manager->getData<std::vector<double>>( ghostID );
    if ( addBufferID != 0 )
        d_AddBuffer = manager->getData<std::vector<double>>( addBufferID );
}


} // namespace AMP::LinearAlgebra


#endif
