#include "AMP/vectors/data/VectorData.h"
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
 * Get/Set values by global id                                   *
 ****************************************************************/
void VectorData::getLocalValuesByGlobalID( size_t N, const size_t *ndx, double *vals ) const
{
    constexpr size_t N_max = 128;
    while ( N > N_max ) {
        getLocalValuesByGlobalID( N_max, ndx, vals );
        N -= N_max;
        ndx  = &ndx[N_max];
        vals = &vals[N_max];
    }
    size_t index[N_max];
    for ( size_t i = 0; i < N; i++ ) {
        AMP_ASSERT( ndx[i] >= d_localStart && ndx[i] < ( d_localStart + d_localSize ) );
        index[i] = ndx[i] - d_localStart;
    }
    getValuesByLocalID( N, index, vals );
}
void VectorData::addLocalValuesByGlobalID( size_t N, const size_t *ndx, const double *vals )
{
    constexpr size_t N_max = 128;
    while ( N > N_max ) {
        addLocalValuesByGlobalID( N_max, ndx, vals );
        N -= N_max;
        ndx  = &ndx[N_max];
        vals = &vals[N_max];
    }
    size_t index[N_max];
    for ( size_t i = 0; i < N; i++ ) {
        AMP_ASSERT( ndx[i] >= d_localStart && ndx[i] < ( d_localStart + d_localSize ) );
        index[i] = ndx[i] - d_localStart;
    }
    addValuesByLocalID( N, index, vals );
}
void VectorData::setLocalValuesByGlobalID( size_t N, const size_t *ndx, const double *vals )
{
    constexpr size_t N_max = 128;
    while ( N > N_max ) {
        setLocalValuesByGlobalID( N_max, ndx, vals );
        N -= N_max;
        ndx  = &ndx[N_max];
        vals = &vals[N_max];
    }
    size_t index[N_max];
    for ( size_t i = 0; i < N; i++ ) {
        AMP_ASSERT( ndx[i] >= d_localStart && ndx[i] < ( d_localStart + d_localSize ) );
        index[i] = ndx[i] - d_localStart;
    }
    setValuesByLocalID( N, index, vals );
}
void VectorData::getValuesByGlobalID( size_t N, const size_t *ndx, double *vals ) const
{
    constexpr size_t N_max = 128;
    while ( N > N_max ) {
        getValuesByGlobalID( N_max, ndx, vals );
        N -= N_max;
        ndx  = &ndx[N_max];
        vals = &vals[N_max];
    }
    size_t N_local = 0, N_ghost = 0;
    size_t local_index[N_max], ghost_index[N_max];
    double local_vals[N_max], ghost_vals[N_max];
    for ( size_t i = 0; i < N; i++ ) {
        if ( ( ndx[i] >= d_localStart ) && ( ndx[i] < ( d_localStart + d_localSize ) ) ) {
            local_index[N_local] = ndx[i] - d_localStart;
            N_local++;
        } else {
            ghost_index[N_ghost] = ndx[i];
            N_ghost++;
        }
    }
    if ( N_local > 0 )
        getValuesByLocalID( N_local, local_index, local_vals );
    if ( N_ghost > 0 )
        getGhostValuesByGlobalID( N_ghost, ghost_index, ghost_vals );
    N_local = 0;
    N_ghost = 0;
    for ( size_t i = 0; i < N; i++ ) {
        if ( ( ndx[i] >= d_localStart ) && ( ndx[i] < ( d_localStart + d_localSize ) ) ) {
            vals[i] = local_vals[N_local];
            N_local++;
        } else {
            vals[i] = ghost_vals[N_ghost];
            N_ghost++;
        }
    }
}
void VectorData::setValuesByGlobalID( size_t N, const size_t *ndx, const double *vals )
{
    constexpr size_t N_max = 128;
    while ( N > N_max ) {
        setValuesByGlobalID( N_max, ndx, vals );
        N -= N_max;
        ndx  = &ndx[N_max];
        vals = &vals[N_max];
    }
    size_t N_local = 0, N_ghost = 0;
    size_t local_index[N_max], ghost_index[N_max];
    double local_vals[N_max], ghost_vals[N_max];
    for ( size_t i = 0; i < N; i++ ) {
        if ( ( ndx[i] >= d_localStart ) && ( ndx[i] < ( d_localStart + d_localSize ) ) ) {
            local_index[N_local] = ndx[i] - d_localStart;
            local_vals[N_local]  = vals[i];
            N_local++;
        } else {
            ghost_index[N_ghost] = ndx[i];
            ghost_vals[N_ghost]  = vals[i];
            N_ghost++;
        }
    }
    if ( N_local > 0 )
        setValuesByLocalID( N_local, local_index, local_vals );
    if ( N_ghost > 0 )
        setGhostValuesByGlobalID( N_ghost, ghost_index, ghost_vals );
}
void VectorData::addValuesByGlobalID( size_t N, const size_t *ndx, const double *vals )
{
    constexpr size_t N_max = 128;
    while ( N > N_max ) {
        addValuesByGlobalID( N_max, ndx, vals );
        N -= N_max;
        ndx  = &ndx[N_max];
        vals = &vals[N_max];
    }
    size_t N_local = 0, N_ghost = 0;
    size_t local_index[N_max], ghost_index[N_max];
    double local_vals[N_max], ghost_vals[N_max];
    for ( size_t i = 0; i < N; i++ ) {
        if ( ( ndx[i] >= d_localStart ) && ( ndx[i] < ( d_localStart + d_localSize ) ) ) {
            local_index[N_local] = ndx[i] - d_localStart;
            local_vals[N_local]  = vals[i];
            N_local++;
        } else {
            ghost_index[N_ghost] = ndx[i];
            ghost_vals[N_ghost]  = vals[i];
            N_ghost++;
        }
    }
    if ( N_local > 0 )
        addValuesByLocalID( N_local, local_index, local_vals );
    if ( N_ghost > 0 )
        addGhostValuesByGlobalID( N_ghost, ghost_index, ghost_vals );
}


/****************************************************************
 * Get/Set ghost values by global id                             *
 ****************************************************************/
void VectorData::setGhostValuesByGlobalID( size_t N, const size_t *ndx, const double *vals )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::ADDING );
    *d_UpdateState = UpdateState::SETTING;
    for ( size_t i = 0; i < N; i++ ) {
        if ( ( ndx[i] < d_localStart ) || ( ndx[i] >= ( d_localStart + d_localSize ) ) ) {
            ( *d_Ghosts )[d_CommList->getLocalGhostID( ndx[i] )] = vals[i];
        } else {
            AMP_ERROR( "Non ghost index" );
        }
    }
}
void VectorData::addGhostValuesByGlobalID( size_t N, const size_t *ndx, const double *vals )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::SETTING );
    *d_UpdateState = UpdateState::ADDING;
    for ( size_t i = 0; i < N; i++ ) {
        if ( ( ndx[i] < d_localStart ) || ( ndx[i] >= ( d_localStart + d_localSize ) ) ) {
            ( *d_AddBuffer )[d_CommList->getLocalGhostID( ndx[i] )] += vals[i];
        } else {
            AMP_ERROR( "Non ghost index" );
        }
    }
}
void VectorData::getGhostValuesByGlobalID( size_t N, const size_t *ndx, double *vals ) const
{
    for ( size_t i = 0; i < N; i++ ) {
        if ( ( ndx[i] < d_localStart ) || ( ndx[i] >= ( d_localStart + d_localSize ) ) ) {
            vals[i] = ( *d_Ghosts )[d_CommList->getLocalGhostID( ndx[i] )] +
                      ( *d_AddBuffer )[d_CommList->getLocalGhostID( ndx[i] )];
        } else {
            AMP_ERROR( "Tried to get a non-ghost ghost value" );
        }
    }
}
void VectorData::getGhostAddValuesByGlobalID( size_t N, const size_t *ndx, double *vals ) const
{
    for ( size_t i = 0; i < N; i++ ) {
        if ( ( ndx[i] < d_localStart ) || ( ndx[i] >= ( d_localStart + d_localSize ) ) ) {
            vals[i] = ( *d_AddBuffer )[d_CommList->getLocalGhostID( ndx[i] )];
        } else {
            AMP_ERROR( "Tried to get a non-ghost ghost value" );
        }
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
            std::vector<double> send_vec_add( d_CommList->getVectorReceiveBufferSize() );
            std::vector<double> recv_vec_add( d_CommList->getVectorSendBufferSize() );
            d_CommList->packReceiveBuffer( send_vec_add, *this );
            d_CommList->scatter_add( send_vec_add, recv_vec_add );
            d_CommList->unpackSendBufferAdd( recv_vec_add, *this );
            for ( auto &elem : *d_AddBuffer ) {
                elem = 0.0;
            }
        }
        *d_UpdateState = UpdateState::SETTING;
        std::vector<double> send_vec( d_CommList->getVectorSendBufferSize() );
        std::vector<double> recv_vec( d_CommList->getVectorReceiveBufferSize() );
        d_CommList->packSendBuffer( send_vec, *this );
        d_CommList->scatter_set( send_vec, recv_vec );
        d_CommList->unpackReceiveBufferSet( recv_vec, *this );
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
 * dump data to ostream                                          *
 ****************************************************************/
void VectorData::copyGhostValues( const VectorData &rhs )
{
    if ( getGhostSize() == 0 ) {
        // No ghosts to fill, we don't need to do anything
    } else if ( getGhostSize() == rhs.getGhostSize() ) {
        // The ghosts in the src vector match the current vector
        // Copy the ghosts from the rhs
        std::vector<size_t> ghostIDs = getCommunicationList()->getGhostIDList();
        std::vector<double> values( ghostIDs.size() );
        rhs.getGhostValuesByGlobalID( ghostIDs.size(), &ghostIDs[0], &values[0] );
        this->setGhostValuesByGlobalID( ghostIDs.size(), &ghostIDs[0], &values[0] );
        // Copy the consistency state from the rhs
        *d_UpdateState = *( rhs.getUpdateStatusPtr() );
    } else {
        // We can't copy the ghosts from the rhs
        // Use makeConsistent to fill the ghosts
        // Note: this will incure global communication
        *d_UpdateState = *( rhs.getUpdateStatusPtr() );
        if ( *d_UpdateState == UpdateState::UNCHANGED )
            *d_UpdateState = UpdateState::LOCAL_CHANGED;
    }
}

/****************************************************************
 * reset a vector                                               *
 ****************************************************************/
void VectorData::reset() { AMP_ERROR( "Not implemented" ); }

void VectorData::print( std::ostream &os, const std::string &name, const std::string &prefix ) const
{
    NULL_USE( os );
    NULL_USE( name );
    NULL_USE( prefix );
    AMP_ERROR( "Not implemented" );
}

} // namespace AMP::LinearAlgebra
