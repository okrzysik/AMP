#include "AMP/vectors/data/VectorData.h"
#include "AMP/IO/RestartManager.h"
#include "AMP/utils/UtilityMacros.h"
#include "AMP/vectors/data/DataChangeListener.h"


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
 * makeConsistent                                                *
 ****************************************************************/
void VectorData::makeConsistent( ScatterType t )
{
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
AMP_MPI VectorData::getComm() const
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

void VectorData::print( std::ostream &os, const std::string &name, const std::string &prefix ) const
{
    NULL_USE( os );
    NULL_USE( name );
    NULL_USE( prefix );
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
void VectorData::registerChildObjects( AMP::IO::RestartManager * ) const
{
    AMP_ERROR( "Need to implement registerChildObjects for " + VectorDataName() );
}
void VectorData::writeRestart( int64_t ) const
{
    AMP_ERROR( "Need to implement writeRestart for " + VectorDataName() );
}


} // namespace AMP::LinearAlgebra


/********************************************************
 *  Restart operations                                   *
 ********************************************************/
template<>
AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::VectorData>::DataStoreType(
    const std::string &name,
    std::shared_ptr<const AMP::LinearAlgebra::VectorData> data,
    RestartManager *manager )
    : d_data( data )
{
    d_name = name;
    d_hash = data->getID();
    d_data->registerChildObjects( manager );
    // Register the communication list
    auto commList = data->getCommunicationList();
    if ( commList )
        manager->registerData( commList );
}
template<>
void AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::VectorData>::write(
    hid_t fid, const std::string &name ) const
{
    hid_t gid = createGroup( fid, name );
    d_data->writeRestart( gid );
    writeHDF5( gid, "ClassType", d_data->VectorDataName() );
    auto commList     = d_data->getCommunicationList();
    auto commListHash = commList ? commList->getID() : 0;
    writeHDF5( gid, "CommListHash", commListHash );
    closeGroup( gid );
}
template<>
std::shared_ptr<AMP::LinearAlgebra::VectorData>
AMP::IO::RestartManager::getData<AMP::LinearAlgebra::VectorData>( const std::string &name )
{
    hid_t gid = openGroup( d_fid, name );
    std::string type;
    readHDF5( gid, "ClassType", type );
    std::shared_ptr<AMP::LinearAlgebra::VectorData> data;
    // Load the object (we will need to replace the if/else with a factory
    AMP_ERROR( "Not finished" );
    /*if ( type == "DOFManager" ) {
        dofs = std::make_shared<AMP::Discretization::DOFManager>( gid, this );
    } else if ( type == "StructuredGeometryMesh" ) {
        // mesh = std::make_shared<AMP::Mesh::StructuredGeometryMesh>( gid, this );
    } else if ( type == "MovableBoxMesh" ) {
        // mesh = std::make_shared<AMP::Mesh::MovableBoxMesh>( gid, this );
    } else {
        AMP_ERROR( "Not finished: " + type );
    }*/
    // Load the communication data
    uint64_t commListHash = 0;
    readHDF5( gid, "CommListHash", commListHash );
    if ( commListHash ) {
        auto commList = getData<AMP::LinearAlgebra::CommunicationList>( commListHash );
        data->setCommunicationList( commList );
    }
    closeGroup( gid );
    return data;
}
