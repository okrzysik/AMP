#include "AMP/vectors/data/VectorData.h"
#include "AMP/IO/RestartManager.h"
#include "AMP/utils/UtilityMacros.h"
#include "AMP/vectors/data/DataChangeListener.h"
#include "AMP/vectors/data/MultiVectorData.h"
#include "AMP/vectors/data/VectorDataDefault.h"
#include "AMP/vectors/data/VectorDataFactory.h"

#include "ProfilerApp.h"


namespace AMP::LinearAlgebra {


VectorData::VectorData()
    : d_UpdateState{ std::make_shared<UpdateState>() },
      d_Ghosts{ std::make_shared<std::vector<double>>() },
      d_AddBuffer{ std::make_shared<std::vector<double>>() }
{
    *d_UpdateState = UpdateState::UNCHANGED;
}

VectorData::VectorData( std::shared_ptr<CommunicationList> comm )
    : d_UpdateState{ std::make_shared<UpdateState>() }
{
    setCommunicationList( comm );
    *d_UpdateState = UpdateState::UNCHANGED;
}

void VectorData::setCommunicationList( std::shared_ptr<CommunicationList> comm )
{
    AMP_ASSERT( comm );
    d_CommList = comm;
    if ( comm ) {
        d_Ghosts =
            std::make_shared<std::vector<double>>( d_CommList->getVectorReceiveBufferSize() );
        d_AddBuffer =
            std::make_shared<std::vector<double>>( d_CommList->getVectorReceiveBufferSize() );
    }
}

/****************************************************************
 * Get/Set ghost values by global id                             *
 ****************************************************************/
void VectorData::setGhostValuesByGlobalID( size_t N,
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
void VectorData::addGhostValuesByGlobalID( size_t N,
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
void VectorData::getGhostValuesByGlobalID( size_t N,
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
void VectorData::getGhostAddValuesByGlobalID( size_t N,
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
 * makeConsistent / UpdateState                                  *
 ****************************************************************/
UpdateState VectorData::getLocalUpdateStatus() const { return *d_UpdateState; }
UpdateState VectorData::getGlobalUpdateStatus() const
{
    PROFILE( "getGlobalUpdateStatus" );
    auto local = getComm().allGather<uint8_t>( static_cast<int8_t>( getLocalUpdateStatus() ) );
    auto state = UpdateState::UNCHANGED;
    for ( auto s : local ) {
        auto s2 = static_cast<UpdateState>( s );
        if ( s2 == UpdateState::UNCHANGED ) {
            continue;
        } else if ( s2 == UpdateState::LOCAL_CHANGED && state == UpdateState::UNCHANGED ) {
            state = UpdateState::LOCAL_CHANGED;
        } else if ( s2 == UpdateState::LOCAL_CHANGED ) {
            continue;
        } else if ( s2 == UpdateState::ADDING &&
                    ( state == UpdateState::UNCHANGED || state == UpdateState::LOCAL_CHANGED ||
                      state == UpdateState::ADDING ) ) {
            state = UpdateState::ADDING;
        } else if ( s2 == UpdateState::SETTING &&
                    ( state == UpdateState::UNCHANGED || state == UpdateState::LOCAL_CHANGED ||
                      state == UpdateState::SETTING ) ) {
            state = UpdateState::SETTING;
        } else {
            state = UpdateState::MIXED;
        }
    }
    return state;
}
void VectorData::setUpdateStatus( UpdateState state ) { *d_UpdateState = state; }
void VectorData::setUpdateStatusPtr( std::shared_ptr<UpdateState> rhs ) { d_UpdateState = rhs; }
std::shared_ptr<UpdateState> VectorData::getUpdateStatusPtr() const { return d_UpdateState; }
void VectorData::makeConsistent( ScatterType t )
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
void VectorData::makeConsistent()
{
    auto state = getGlobalUpdateStatus();
    if ( state == UpdateState::UNCHANGED ) {
        return;
    } else if ( state == UpdateState::LOCAL_CHANGED || state == UpdateState::SETTING ) {
        makeConsistent( ScatterType::CONSISTENT_SET );
    } else if ( state == UpdateState::ADDING || state == UpdateState::MIXED ) {
        makeConsistent( ScatterType::CONSISTENT_ADD );
    } else {
        AMP_ERROR( "Unknown update status" );
    }
}


/****************************************************************
 * dataChanged                                                   *
 ****************************************************************/
void VectorData::dataChanged()
{
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
    fireDataChange();
}


/****************************************************************
 * Default clone                                                 *
 ****************************************************************/
const AMP_MPI &VectorData::getComm() const
{
    AMP_ASSERT( d_CommList );
    return d_CommList->getComm();
}


/****************************************************************
 * Check if two data blocks are alias to each other              *
 ****************************************************************/
bool VectorData::isAnAliasOf( const VectorData &rhs ) const
{
    if ( numberOfDataBlocks() != rhs.numberOfDataBlocks() )
        return false;
    for ( size_t i = 0; i < numberOfDataBlocks(); i++ ) {
        if ( sizeOfDataBlock( i ) != rhs.sizeOfDataBlock( i ) ||
             getRawDataBlockAsVoid( i ) != rhs.getRawDataBlockAsVoid( i ) )
            return false;
    }
    return true;
}


/****************************************************************
 * dump data to ostream                                          *
 ****************************************************************/
void VectorData::dumpOwnedData( std::ostream &out, size_t GIDoffset, size_t LIDoffset ) const
{
    auto curElement = begin();
    size_t gid      = GIDoffset;
    if ( getCommunicationList() )
        gid += getCommunicationList()->getStartGID();
    size_t lid = LIDoffset;
    while ( curElement != end() ) {
        out << "  GID: " << gid << "  LID: " << lid << "  Value: " << *curElement << "\n";
        ++curElement;
        ++gid;
        ++lid;
    }
}
void VectorData::dumpGhostedData( std::ostream &out, size_t offset ) const
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
 * Component data                                                *
 ****************************************************************/
size_t VectorData::getNumberOfComponents() const { return 1; }
std::shared_ptr<VectorData> VectorData::getComponent( size_t i )
{
    AMP_ASSERT( getNumberOfComponents() == 1 );
    AMP_ASSERT( i == 0 );
    return shared_from_this();
}
std::shared_ptr<const VectorData> VectorData::getComponent( size_t i ) const
{
    AMP_ASSERT( getNumberOfComponents() == 1 );
    AMP_ASSERT( i == 0 );
    return shared_from_this();
}


/****************************************************************
 * dump data to ostream                                          *
 ****************************************************************/
void VectorData::copyGhostValues( const VectorData &rhs )
{
    if ( getGhostSize() == 0 ) {
        // No ghosts to fill, copy the consistency state from the rhs
        *d_UpdateState = *rhs.getUpdateStatusPtr();
    } else if ( getGhostSize() == rhs.getGhostSize() ) {
        // The ghosts in the src vector match the current vector
        // Copy the ghosts from the rhs
        auto ghostIDs = getCommunicationList()->getGhostIDList();
        std::vector<double> values( ghostIDs.size() );
        rhs.getGhostValuesByGlobalID( ghostIDs.size(), &ghostIDs[0], &values[0] );
        this->setGhostValuesByGlobalID( ghostIDs.size(), &ghostIDs[0], &values[0] );
        // Copy the consistency state from the rhs
        *d_UpdateState = *rhs.getUpdateStatusPtr();
    } else {
        // We can't copy the ghosts from the rhs
        // Use makeConsistent to fill the ghosts
        // Note: this will insure global communication
        *d_UpdateState = *rhs.getUpdateStatusPtr();
        if ( *d_UpdateState == UpdateState::UNCHANGED )
            *d_UpdateState = UpdateState::LOCAL_CHANGED;
    }
}

/****************************************************************
 * reset vector data                                             *
 ****************************************************************/
void VectorData::reset() {}

void VectorData::print( std::ostream &, const std::string &, const std::string & ) const
{
    AMP_ERROR( "Not implemented" );
}


/****************************************************************
 * Get an id                                                     *
 ****************************************************************/
uint64_t VectorData::getID() const
{
    return getComm().bcast( reinterpret_cast<uint64_t>( this ), 0 );
}


/****************************************************************
 * Write/Read restart data                                       *
 ****************************************************************/
void VectorData::registerChildObjects( AMP::IO::RestartManager *manager ) const
{
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
void VectorData::writeRestart( int64_t fid ) const
{
    uint64_t commListID  = d_CommList ? d_CommList->getID() : 0;
    uint64_t updateID    = reinterpret_cast<uint64_t>( d_UpdateState.get() );
    uint64_t ghostID     = reinterpret_cast<uint64_t>( d_Ghosts.get() );
    uint64_t addBufferID = reinterpret_cast<uint64_t>( d_AddBuffer.get() );
    IO::writeHDF5( fid, "localSize", d_localSize );
    IO::writeHDF5( fid, "globalSize", d_globalSize );
    IO::writeHDF5( fid, "localStart", d_localStart );
    IO::writeHDF5( fid, "commListID", commListID );
    IO::writeHDF5( fid, "updateID", updateID );
    IO::writeHDF5( fid, "ghostID", ghostID );
    IO::writeHDF5( fid, "addBufferID", addBufferID );
}


VectorData::VectorData( int64_t fid, AMP::IO::RestartManager *manager )
{
    uint64_t commListID, updateID, ghostID, addBufferID;
    IO::readHDF5( fid, "localSize", d_localSize );
    IO::readHDF5( fid, "globalSize", d_globalSize );
    IO::readHDF5( fid, "localStart", d_localStart );
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


/********************************************************
 *  Restart operations                                   *
 ********************************************************/
template<>
AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::VectorData>::DataStoreType(
    std::shared_ptr<const AMP::LinearAlgebra::VectorData> data, RestartManager *manager )
    : d_data( data )
{
    d_hash = data->getID();
    d_data->registerChildObjects( manager );
    // Register the communication list
    auto commList = data->getCommunicationList();
    if ( commList )
        manager->registerObject( commList );
}
template<>
void AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::VectorData>::write(
    hid_t fid, const std::string &name ) const
{
    hid_t gid         = createGroup( fid, name );
    auto commList     = d_data->getCommunicationList();
    auto commListHash = commList ? commList->getID() : 0;
    writeHDF5( gid, "ClassType", d_data->VectorDataName() );
    writeHDF5( gid, "CommListHash", commListHash );
    d_data->writeRestart( gid );
    closeGroup( gid );
}
template<>
std::shared_ptr<AMP::LinearAlgebra::VectorData>
AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::VectorData>::read(
    hid_t fid, const std::string &name, RestartManager *manager ) const
{
    hid_t gid    = openGroup( fid, name );
    auto vecData = AMP::LinearAlgebra::VectorDataFactory::create( gid, manager );
    closeGroup( gid );
    return vecData;
}
template<>
AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::UpdateState>::DataStoreType(
    std::shared_ptr<const AMP::LinearAlgebra::UpdateState> data, RestartManager * )
    : d_data( data )
{
    d_hash = reinterpret_cast<uint64_t>( d_data.get() );
}
template<>
void AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::UpdateState>::write(
    hid_t fid, const std::string &name ) const
{
    int value = static_cast<int>( *d_data );
    writeHDF5( fid, name, value );
}
template<>
std::shared_ptr<AMP::LinearAlgebra::UpdateState>
AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::UpdateState>::read(
    hid_t fid, const std::string &name, RestartManager * ) const
{
    int value;
    readHDF5( fid, name, value );
    auto state = static_cast<AMP::LinearAlgebra::UpdateState>( value );
    return std::make_shared<AMP::LinearAlgebra::UpdateState>( state );
}
