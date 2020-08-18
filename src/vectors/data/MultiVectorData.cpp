#ifndef included_AMP_MultiVectorData_hpp
#define included_AMP_MultiVectorData_hpp

#include "AMP/vectors/data/MultiVectorData.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/vectors/data/VectorData.h"
#include "AMP/vectors/Vector.h"

#include "ProfilerApp.h"


namespace AMP {
namespace LinearAlgebra {


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
size_t MultiVectorData::getLocalSize() const { return d_globalDOFManager->numLocalDOF(); }
size_t MultiVectorData::getGlobalSize() const { return d_globalDOFManager->numGlobalDOF(); }
size_t MultiVectorData::getGhostSize() const
{
    size_t ans = 0;
    for ( const auto &data : d_data )
        ans += data->getGhostSize();
    return ans;
}

uint64_t MultiVectorData::getDataID() const { return 0; }
bool MultiVectorData::isTypeId( size_t hash, size_t block ) const
{
    size_t curOffset = 0;
    for ( const auto &data : d_data ) {
        curOffset += data->numberOfDataBlocks();
        if ( block < curOffset ) {
            size_t index = block + data->numberOfDataBlocks() - curOffset;
            return data->isTypeId( hash, index );
        }
    }
    return false;
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
void MultiVectorData::setValuesByLocalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionLocalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->setValuesByLocalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVectorData::setLocalValuesByGlobalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->setLocalValuesByGlobalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVectorData::setGhostValuesByGlobalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->setGhostValuesByGlobalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVectorData::setValuesByGlobalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->setValuesByGlobalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVectorData::addValuesByLocalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionLocalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->addValuesByLocalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVectorData::addLocalValuesByGlobalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->addLocalValuesByGlobalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVectorData::addValuesByGlobalID( int num, size_t *indices, const double *in_vals )
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    partitionGlobalValues( num, indices, in_vals, ndxs, vals );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->addValuesByGlobalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
}
void MultiVectorData::getValuesByGlobalID( int num, size_t *indices, double *out_vals ) const
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    std::vector<std::vector<int>> remap;
    partitionGlobalValues( num, indices, out_vals, ndxs, vals, &remap );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->getValuesByGlobalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
    for ( size_t i = 0; i != remap.size(); i++ ) {
        for ( size_t j = 0; j != remap[i].size(); j++ )
            out_vals[remap[i][j]] = vals[i][j];
    }
}
void MultiVectorData::getLocalValuesByGlobalID( int num, size_t *indices, double *out_vals ) const
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    std::vector<std::vector<int>> remap;
    partitionGlobalValues( num, indices, out_vals, ndxs, vals, &remap );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->getLocalValuesByGlobalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
    for ( size_t i = 0; i != remap.size(); i++ ) {
        for ( size_t j = 0; j != remap[i].size(); j++ )
            out_vals[remap[i][j]] = vals[i][j];
    }
}
void MultiVectorData::getGhostValuesByGlobalID( int num, size_t *indices, double *out_vals ) const
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    std::vector<std::vector<int>> remap;
    partitionGlobalValues( num, indices, out_vals, ndxs, vals, &remap );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->getGhostValuesByGlobalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
    for ( size_t i = 0; i != remap.size(); i++ ) {
        for ( size_t j = 0; j != remap[i].size(); j++ )
            out_vals[remap[i][j]] = vals[i][j];
    }
}
void MultiVectorData::getValuesByLocalID( int num, size_t *indices, double *out_vals ) const
{
    if ( num == 0 )
        return;
    std::vector<std::vector<size_t>> ndxs;
    std::vector<std::vector<double>> vals;
    std::vector<std::vector<int>> remap;
    partitionLocalValues( num, indices, out_vals, ndxs, vals, &remap );
    for ( size_t i = 0; i != ndxs.size(); i++ ) {
        if ( ndxs[i].size() )
            d_data[i]->getValuesByLocalID( ndxs[i].size(), &( ndxs[i][0] ), &( vals[i][0] ) );
    }
    for ( size_t i = 0; i != remap.size(); i++ ) {
        for ( size_t j = 0; j != remap[i].size(); j++ )
            out_vals[remap[i][j]] = vals[i][j];
    }
}


/****************************************************************
 * Copy raw data                                                 *
 ****************************************************************/
void MultiVectorData::putRawData( const double *in )
{
    int cur_off = 0;
    for ( const auto &data : d_data ) {
        data->putRawData( in + cur_off );
        cur_off += data->getLocalSize();
    }
}

void MultiVectorData::copyOutRawData( double *out ) const
{
    size_t curOffset = 0;
    for ( const auto &data : d_data ) {
        data->copyOutRawData( &out[curOffset] );
        curOffset += data->getLocalSize();
    }
}


/****************************************************************
 * makeConsistent                                                *
 ****************************************************************/
void MultiVectorData::makeConsistent( ScatterType t )
{
    for ( const auto &data : d_data )
        data->makeConsistent( t );
    *d_UpdateState = VectorData::UpdateState::UNCHANGED;
}
VectorData::UpdateState MultiVectorData::getUpdateStatus() const
{
    VectorData::UpdateState state = *d_UpdateState;
    for ( const auto &data : d_data ) {
        VectorData::UpdateState sub_state = data->getUpdateStatus();
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
 * Function to partition the global ids by the sub vectors       *
 ****************************************************************/
void MultiVectorData::partitionGlobalValues( const int num,
                                             const size_t *indices,
                                             const double *vals,
                                             std::vector<std::vector<size_t>> &out_indices,
                                             std::vector<std::vector<double>> &out_vals,
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
    auto *manager = dynamic_cast<AMP::Discretization::multiDOFManager *>( d_globalDOFManager );
    for ( size_t i = 0; i < d_data.size(); i++ ) {
        std::vector<size_t> subDOFs = manager->getSubDOF( i, globalDOFs );
        size_t count                = 0;
        for ( auto &subDOF : subDOFs ) {
            if ( subDOF != neg_one )
                count++;
        }
        out_indices[i] = std::vector<size_t>( count, neg_one );
        out_vals[i]    = std::vector<double>( count, 0.0 );
        if ( remap != nullptr )
            remap->operator[]( i ) = std::vector<int>( count, -1 );
        count = 0;
        for ( size_t j = 0; j < subDOFs.size(); j++ ) {
            if ( subDOFs[j] != neg_one ) {
                out_indices[i][count] = subDOFs[j];
                out_vals[i][count]    = vals[j];
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
                                            const double *vals,
                                            std::vector<std::vector<size_t>> &out_indices,
                                            std::vector<std::vector<double>> &out_vals,
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
    partitionGlobalValues( num, &global_indices[0], vals, out_indices, out_vals, remap );
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
  auto vec = dynamic_cast<Vector *>(d_data[i]);
  if (vec) return vec->getVectorData();
  return d_data[i];
}

const VectorData *MultiVectorData::getVectorData( size_t i ) const
{
  auto vec = dynamic_cast<Vector const *>(d_data[i]);
  if (vec) return vec->getVectorData();
  return d_data[i];
}

} // namespace LinearAlgebra
} // namespace AMP

#endif
