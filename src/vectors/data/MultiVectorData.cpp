#ifndef included_AMP_MultiVectorData_hpp
#define included_AMP_MultiVectorData_hpp

#include "AMP/vectors/data/MultiVectorData.h"
#include "AMP/IO/RestartManager.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorData.h"

#include "ProfilerApp.h"

#include <cstddef>


namespace AMP::LinearAlgebra {

void MultiVectorData::resetMultiVectorData( AMP::Discretization::DOFManager *manager,
                                            const std::vector<VectorData *> &data )
{
    d_data = data;
    AMP_ASSERT( manager );
    d_globalDOFManager = manager;
    auto globalMgr     = dynamic_cast<AMP::Discretization::multiDOFManager *>( d_globalDOFManager );
    AMP_ASSERT( globalMgr );
    auto subManagers = globalMgr->getDOFManagers();
    d_subDOFManager.resize( subManagers.size() );
    for ( size_t i = 0; i < subManagers.size(); i++ )
        d_subDOFManager[i] = subManagers[i].get();

    // Create a new communication list
    auto remote_DOFs = globalMgr->getRemoteDOFs();
    bool ghosts      = globalMgr->getComm().anyReduce( !remote_DOFs.empty() );
    if ( !ghosts ) {
        d_CommList =
            std::make_shared<CommunicationList>( globalMgr->numLocalDOF(), globalMgr->getComm() );
    } else {
        auto params           = std::make_shared<AMP::LinearAlgebra::CommunicationListParameters>();
        params->d_comm        = globalMgr->getComm();
        params->d_localsize   = globalMgr->numLocalDOF();
        params->d_remote_DOFs = remote_DOFs;
        d_CommList            = std::make_shared<AMP::LinearAlgebra::CommunicationList>( params );
    }

    // Initialize local/global size
    d_localSize  = d_globalDOFManager->numLocalDOF();
    d_globalSize = d_globalDOFManager->numGlobalDOF();
    d_localStart = d_CommList->getStartGID();
}

void MultiVectorData::reset()
{
    for ( auto &vdata : d_data )
        vdata->reset();
}

/****************************************************************
 * Return basic properties                                       *
 ****************************************************************/
size_t MultiVectorData::numberOfDataBlocks() const
{
    size_t ans = 0;
    for ( const auto &data : d_data )
        ans += data->numberOfDataBlocks();
    return ans;
}
size_t MultiVectorData::sizeOfDataBlock( size_t i ) const
{
    size_t retVal = 0;
    size_t rightOffset, leftOffset;
    rightOffset = leftOffset = 0;
    for ( const auto &data : d_data ) {
        rightOffset += data->numberOfDataBlocks();
        if ( i < rightOffset ) {
            retVal = data->sizeOfDataBlock( i - leftOffset );
            break;
        }
        leftOffset = rightOffset;
    }
    return retVal;
}
size_t MultiVectorData::getGhostSize() const
{
    size_t ans = 0;
    for ( const auto &data : d_data )
        ans += data->getGhostSize();
    return ans;
}
uint64_t MultiVectorData::getDataID() const { return 0; }
typeID MultiVectorData::getType( size_t block ) const
{
    size_t curOffset = 0;
    for ( const auto &data : d_data ) {
        curOffset += data->numberOfDataBlocks();
        if ( block < curOffset ) {
            size_t index = block + data->numberOfDataBlocks() - curOffset;
            return data->getType( index );
        }
    }
    return typeID();
}
size_t MultiVectorData::sizeofDataBlockType( size_t block ) const
{
    size_t curOffset = 0;
    for ( const auto &data : d_data ) {
        curOffset += data->numberOfDataBlocks();
        if ( block < curOffset ) {
            size_t index = block + data->numberOfDataBlocks() - curOffset;
            return data->sizeofDataBlockType( index );
        }
    }
    return 0;
}


/****************************************************************
 * Component data                                                *
 ****************************************************************/
size_t MultiVectorData::getNumberOfComponents() const
{
    size_t N = 0;
    for ( auto data : d_data )
        N += data->getNumberOfComponents();
    return N;
}
std::shared_ptr<VectorData> MultiVectorData::getComponent( size_t i )
{
    for ( auto data : d_data ) {
        size_t N = data->getNumberOfComponents();
        if ( i < N )
            return data->getComponent( i );
        i -= N;
    }
    return nullptr;
}
std::shared_ptr<const VectorData> MultiVectorData::getComponent( size_t i ) const
{
    for ( auto data : d_data ) {
        size_t N = data->getNumberOfComponents();
        if ( i < N )
            return data->getComponent( i );
        i -= N;
    }
    return nullptr;
}


/****************************************************************
 * Access the raw data blocks                                    *
 ****************************************************************/
void *MultiVectorData::getRawDataBlockAsVoid( size_t i )
{
    size_t curOffset = 0;
    for ( const auto &data : d_data ) {
        curOffset += data->numberOfDataBlocks();
        if ( i < curOffset ) {
            size_t index = i + data->numberOfDataBlocks() - curOffset;
            return data->getRawDataBlock<double>( index );
        }
    }
    return nullptr;
}
const void *MultiVectorData::getRawDataBlockAsVoid( size_t i ) const
{
    size_t curOffset = 0;
    for ( const auto &data : d_data ) {
        curOffset += data->numberOfDataBlocks();
        if ( i < curOffset ) {
            size_t index = i + data->numberOfDataBlocks() - curOffset;
            return data->getRawDataBlock<double>( index );
        }
    }
    return nullptr;
}


/****************************************************************
 * Functions to access data by ID                                *
 ****************************************************************/
void MultiVectorData::setValuesByLocalID( size_t N,
                                          const size_t *indices,
                                          const void *in_vals,
                                          const typeID &id )
{
    if ( N == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<std::byte>> vals;
    partitionLocalValues( N, indices, in_vals, id.bytes, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->setValuesByLocalID( ndxs[i].size(), ndxs[i].data(), vals[i].data(), id );
    }
}
void MultiVectorData::addValuesByLocalID( size_t N,
                                          const size_t *indices,
                                          const void *in_vals,
                                          const typeID &id )
{
    if ( N == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<std::byte>> vals;
    partitionLocalValues( N, indices, in_vals, id.bytes, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->addValuesByLocalID( ndxs[i].size(), ndxs[i].data(), vals[i].data(), id );
    }
}
void MultiVectorData::getValuesByLocalID( size_t N,
                                          const size_t *indices,
                                          void *out_vals,
                                          const typeID &id ) const
{
    if ( N == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<std::byte>> vals;
    std::vector<std::vector<int>> remap;
    partitionLocalValues( N, indices, out_vals, id.bytes, ndxs, vals, &remap );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->getValuesByLocalID( ndxs[i].size(), ndxs[i].data(), vals[i].data(), id );
    }
    auto out = reinterpret_cast<std::byte *>( out_vals );
    for ( size_t i = 0; i != remap.size(); i++ ) {
        for ( size_t j = 0; j != remap[i].size(); j++ )
            memcpy( &out[remap[i][j] * id.bytes], &vals[i][j * id.bytes], id.bytes );
    }
}
void MultiVectorData::setGhostValuesByGlobalID( size_t N,
                                                const size_t *indices,
                                                const void *in_vals,
                                                const typeID &id )
{
    if ( N == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<std::byte>> vals;
    partitionGlobalValues( N, indices, in_vals, id.bytes, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->setGhostValuesByGlobalID(
                ndxs[i].size(), ndxs[i].data(), vals[i].data(), id );
    }
}
void MultiVectorData::addGhostValuesByGlobalID( size_t N,
                                                const size_t *indices,
                                                const void *in_vals,
                                                const typeID &id )
{
    if ( N == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<std::byte>> vals;
    partitionGlobalValues( N, indices, in_vals, id.bytes, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->addGhostValuesByGlobalID(
                ndxs[i].size(), ndxs[i].data(), vals[i].data(), id );
    }
}
void MultiVectorData::getGhostValuesByGlobalID( size_t N,
                                                const size_t *indices,
                                                void *out_vals,
                                                const typeID &id ) const
{
    if ( N == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<std::byte>> vals;
    std::vector<std::vector<int>> remap;
    partitionGlobalValues( N, indices, out_vals, id.bytes, ndxs, vals, &remap );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->getGhostValuesByGlobalID(
                ndxs[i].size(), ndxs[i].data(), vals[i].data(), id );
    }
    auto out = reinterpret_cast<std::byte *>( out_vals );
    for ( size_t i = 0; i != remap.size(); i++ ) {
        for ( size_t j = 0; j != remap[i].size(); j++ )
            memcpy( &out[remap[i][j] * id.bytes], &vals[i][j * id.bytes], id.bytes );
    }
}


/****************************************************************
 * Copy raw data                                                 *
 ****************************************************************/
void MultiVectorData::putRawData( const void *in, const typeID &id )
{
    const char *ptr = reinterpret_cast<const char *>( in );
    for ( const auto &data : d_data ) {
        data->putRawData( ptr, id );
        ptr += id.bytes * data->getLocalSize();
    }
}

void MultiVectorData::getRawData( void *out, const typeID &id ) const
{
    char *ptr = reinterpret_cast<char *>( out );
    for ( const auto &data : d_data ) {
        data->getRawData( ptr, id );
        ptr += id.bytes * data->getLocalSize();
    }
}


/****************************************************************
 * makeConsistent                                                *
 ****************************************************************/
void MultiVectorData::makeConsistent( ScatterType t )
{
    for ( const auto &data : d_data ) {
        auto vec = dynamic_cast<Vector *>( data );
        if ( vec ) {
            vec->getVectorData()->makeConsistent( t );
        } else {
            data->makeConsistent( t );
        }
    }
    *d_UpdateState = UpdateState::UNCHANGED;
}
UpdateState MultiVectorData::getUpdateStatus() const
{
    UpdateState state = *d_UpdateState;
    for ( const auto &data : d_data ) {
        UpdateState sub_state = data->getUpdateStatus();
        if ( sub_state == UpdateState::UNCHANGED ) {
            continue;
        } else if ( sub_state == UpdateState::LOCAL_CHANGED && state == UpdateState::UNCHANGED ) {
            state = UpdateState::LOCAL_CHANGED;
        } else if ( sub_state == UpdateState::LOCAL_CHANGED ) {
            continue;
        } else if ( sub_state == UpdateState::ADDING &&
                    ( state == UpdateState::UNCHANGED || state == UpdateState::LOCAL_CHANGED ||
                      state == UpdateState::ADDING ) ) {
            state = UpdateState::ADDING;
        } else if ( sub_state == UpdateState::SETTING &&
                    ( state == UpdateState::UNCHANGED || state == UpdateState::LOCAL_CHANGED ||
                      state == UpdateState::SETTING ) ) {
            state = UpdateState::SETTING;
        } else {
            state = UpdateState::MIXED;
        }
    }
    return state;
}
void MultiVectorData::setUpdateStatus( UpdateState state )
{
    *d_UpdateState = state;
    for ( const auto &data : d_data )
        data->setUpdateStatus( state );
}


/****************************************************************
 * Swap raw data                                                 *
 ****************************************************************/
void MultiVectorData::swapData( VectorData & ) { AMP_ERROR( "Not finished" ); }


/****************************************************************
 * Clone raw data                                                *
 ****************************************************************/
std::shared_ptr<VectorData> MultiVectorData::cloneData( const std::string & ) const
{
    AMP_ERROR( "Not finished" );
    return std::shared_ptr<VectorData>();
}


/****************************************************************
 * Function to partition the global ids by the sub vectors       *
 ****************************************************************/
void MultiVectorData::partitionGlobalValues( const int num,
                                             const size_t *indices,
                                             const void *vals,
                                             size_t bytes,
                                             std::vector<std::vector<size_t>> &out_indices,
                                             std::vector<std::vector<std::byte>> &out_vals,
                                             std::vector<std::vector<int>> *remap ) const
{
    PROFILE_START( "partitionGlobalValues", 2 );
    const size_t neg_one = ~( (size_t) 0 );
    std::vector<size_t> globalDOFs( num, neg_one );
    for ( int i = 0; i < num; i++ )
        globalDOFs[i] = indices[i];
    out_indices.resize( d_data.size() );
    out_vals.resize( d_data.size() );
    if ( remap != nullptr )
        remap->resize( d_data.size() );
    auto data     = reinterpret_cast<const std::byte *>( vals );
    auto *manager = dynamic_cast<AMP::Discretization::multiDOFManager *>( d_globalDOFManager );
    for ( size_t i = 0; i < d_data.size(); i++ ) {
        auto subDOFs = manager->getSubDOF( i, globalDOFs );
        size_t count = 0;
        for ( auto &subDOF : subDOFs ) {
            if ( subDOF != neg_one )
                count++;
        }
        out_indices[i].resize( count, neg_one );
        out_vals[i].resize( count * bytes );
        if ( remap != nullptr )
            remap->operator[]( i ) = std::vector<int>( count, -1 );
        count = 0;
        for ( size_t j = 0; j < subDOFs.size(); j++ ) {
            if ( subDOFs[j] != neg_one ) {
                out_indices[i][count] = subDOFs[j];
                memcpy( &out_vals[i][count * bytes], &data[j * bytes], bytes );
                if ( remap != nullptr )
                    remap->operator[]( i )[count] = j;
                count++;
            }
        }
    }
    PROFILE_STOP( "partitionGlobalValues", 2 );
}


/****************************************************************
 * Function to partition the local ids by the sub vectors       *
 ****************************************************************/
void MultiVectorData::partitionLocalValues( const int num,
                                            const size_t *indices,
                                            const void *vals,
                                            size_t bytes,
                                            std::vector<std::vector<size_t>> &out_indices,
                                            std::vector<std::vector<std::byte>> &out_vals,
                                            std::vector<std::vector<int>> *remap ) const
{
    if ( num == 0 )
        return;
    PROFILE_START( "partitionLocalValues", 2 );
    // Convert the local ids to global ids
    size_t begin_DOF = d_globalDOFManager->beginDOF();
    size_t end_DOF   = d_globalDOFManager->endDOF();
    std::vector<size_t> global_indices( num );
    for ( int i = 0; i < num; i++ ) {
        AMP_INSIST( indices[i] < end_DOF, "Invalid local id" );
        global_indices[i] = indices[i] + begin_DOF;
    }
    // Partition based on the global ids
    partitionGlobalValues( num, &global_indices[0], vals, bytes, out_indices, out_vals, remap );
    // Convert the new global ids back to local ids
    const size_t neg_one = ~( (size_t) 0 );
    for ( size_t i = 0; i < d_data.size(); i++ ) {
        if ( out_indices[i].size() == 0 )
            continue;
        begin_DOF = d_subDOFManager[i]->beginDOF();
        end_DOF   = d_subDOFManager[i]->endDOF();
        for ( auto &elem : out_indices[i] ) {
            AMP_ASSERT( elem != neg_one );
            elem -= begin_DOF;
            AMP_ASSERT( elem < end_DOF );
        }
    }
    PROFILE_STOP( "partitionLocalValues", 2 );
}

VectorData *MultiVectorData::getVectorData( size_t i )
{
    auto vec = dynamic_cast<Vector *>( d_data[i] );
    if ( vec )
        return vec->getVectorData().get();
    return d_data[i];
}

const VectorData *MultiVectorData::getVectorData( size_t i ) const
{
    auto vec = dynamic_cast<Vector const *>( d_data[i] );
    if ( vec )
        return vec->getVectorData().get();
    return d_data[i];
}

/****************************************************************
 * Functions to print the data                                   *
 ****************************************************************/
void MultiVectorData::assemble()
{
    for ( auto vec : d_data )
        vec->assemble();
}


/****************************************************************
 * Functions to print the data                                   *
 ****************************************************************/
void MultiVectorData::dumpOwnedData( std::ostream &out, size_t GIDoffset, size_t LIDoffset ) const
{
    size_t localOffset = 0;
    auto *manager      = (AMP::Discretization::multiDOFManager *) d_globalDOFManager;
    AMP_ASSERT( manager->getDOFManagers().size() == d_subDOFManager.size() );
    for ( size_t i = 0; i != d_subDOFManager.size(); i++ ) {
        auto subManager = d_subDOFManager[i];
        std::vector<size_t> subStartDOF( 1, subManager->beginDOF() );
        auto globalStartDOF = manager->getGlobalDOF( i, subStartDOF );
        size_t globalOffset = globalStartDOF[0] - subStartDOF[0];
        d_data[i]->dumpOwnedData( out, GIDoffset + globalOffset, LIDoffset + localOffset );
        localOffset += d_data[i]->getLocalSize();
    }
}
void MultiVectorData::dumpGhostedData( std::ostream &out, size_t offset ) const
{
    auto manager = (AMP::Discretization::multiDOFManager *) d_globalDOFManager;
    AMP_ASSERT( manager->getDOFManagers().size() == d_subDOFManager.size() );
    for ( size_t i = 0; i != d_data.size(); i++ ) {
        auto subManager = d_subDOFManager[i];
        std::vector<size_t> subStartDOF( 1, subManager->beginDOF() );
        auto globalStartDOF = manager->getGlobalDOF( i, subStartDOF );
        size_t globalOffset = globalStartDOF[0] - subStartDOF[0];
        d_data[i]->dumpGhostedData( out, offset + globalOffset );
    }
}


/****************************************************************
 * Write/Read restart data                                       *
 ****************************************************************/
void MultiVectorData::registerChildObjects( AMP::IO::RestartManager *manager ) const
{
    manager->registerComm( d_comm );
    for ( auto data : d_data )
        manager->registerObject( data->shared_from_this() );
    manager->registerObject( d_globalDOFManager->shared_from_this() );
    for ( auto dofs : d_subDOFManager )
        manager->registerObject( dofs->shared_from_this() );
}
void MultiVectorData::writeRestart( int64_t fid ) const
{
    std::vector<uint64_t> dataHash( d_data.size() );
    std::vector<uint64_t> dofsHash( d_subDOFManager.size() );
    for ( size_t i = 0; i < d_data.size(); i++ )
        dataHash[i] = d_data[i]->getID();
    for ( size_t i = 0; i < d_subDOFManager.size(); i++ )
        dofsHash[i] = d_subDOFManager[i]->getID();
    writeHDF5( fid, "CommHash", d_comm.hash() );
    writeHDF5( fid, "globalDOFsHash", d_globalDOFManager->getID() );
    writeHDF5( fid, "VectorDataHash", dataHash );
    writeHDF5( fid, "DOFManagerHash", dofsHash );
}
MultiVectorData::MultiVectorData( int64_t fid, AMP::IO::RestartManager *manager )
{
    uint64_t commHash, globalDOFsHash;
    std::vector<uint64_t> vectorDataHash, DOFManagerHash;
    readHDF5( fid, "CommHash", commHash );
    readHDF5( fid, "globalDOFsHash", globalDOFsHash );
    readHDF5( fid, "VectorDataHash", vectorDataHash );
    readHDF5( fid, "DOFManagerHash", DOFManagerHash );
    d_comm = manager->getComm( commHash );
    d_data.resize( vectorDataHash.size() );
    for ( size_t i = 0; i < d_data.size(); i++ )
        d_data[i] = manager->getData<VectorData>( vectorDataHash[i] ).get();
    d_globalDOFManager = manager->getData<AMP::Discretization::DOFManager>( globalDOFsHash ).get();
    d_subDOFManager.resize( DOFManagerHash.size() );
    for ( size_t i = 0; i < d_subDOFManager.size(); i++ )
        d_subDOFManager[i] =
            manager->getData<AMP::Discretization::DOFManager>( DOFManagerHash[i] ).get();
}


} // namespace AMP::LinearAlgebra

#endif
