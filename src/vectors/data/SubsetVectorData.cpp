#include "AMP/vectors/data/SubsetVectorData.h"
#include "AMP/discretization/subsetDOFManager.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"

#include "ProfilerApp.h"

#include <algorithm>


namespace AMP::LinearAlgebra {


/****************************************************************
 * Constructors                                                   *
 ****************************************************************/
SubsetVectorData::SubsetVectorData( std::shared_ptr<SubsetVectorParameters> params )
    : VectorData( params->d_CommList ),
      d_ViewVector( params->d_ViewVector ),
      d_DOFManager{ params->d_DOFManager }
{
    auto parentDOF       = d_ViewVector->getDOFManager();
    d_parentLocalStartID = parentDOF->beginDOF();
    auto tmp = std::dynamic_pointer_cast<AMP::Discretization::subsetDOFManager>( d_DOFManager );
    if ( tmp != nullptr ) {
        d_SubsetLocalIDToViewGlobalID = tmp->getLocalParentDOFs();
    } else if ( d_DOFManager->numLocalDOF() == parentDOF->numLocalDOF() ) {
        d_SubsetLocalIDToViewGlobalID.resize( parentDOF->numLocalDOF() );
        size_t beginDOF = parentDOF->beginDOF();
        size_t endDOF   = parentDOF->endDOF();
        for ( size_t i = beginDOF; i < endDOF; i++ )
            d_SubsetLocalIDToViewGlobalID[i - beginDOF] = i;
    } else {
        AMP_ERROR( "Internal error with SubsetVector" );
    }
    // Get a pointer to every value in the subset
    std::vector<double *> data_ptr( d_SubsetLocalIDToViewGlobalID.size(), nullptr );
    auto iterator   = d_ViewVector->constBegin();
    size_t last_pos = d_ViewVector->getCommunicationList()->getStartGID();
    for ( size_t i = 0; i < data_ptr.size(); i++ ) {
        iterator += (int) ( d_SubsetLocalIDToViewGlobalID[i] - last_pos );
        last_pos    = d_SubsetLocalIDToViewGlobalID[i];
        data_ptr[i] = const_cast<double *>( &( *iterator ) );
    }
    // Create the data blocks
    // For now use one datablock for each value, this needs to be changed
    d_dataBlockPtr  = data_ptr;
    d_dataBlockSize = std::vector<size_t>( data_ptr.size(), 1 );
    // Cache the local/global size
    AMP_ASSERT( d_ViewVector );
    d_localSize  = d_CommList->numLocalRows();
    d_globalSize = d_CommList->getTotalSize();
    d_localStart = d_CommList->getStartGID();
}

/****************************************************************
 * Functions to access the raw data blocks                       *
 ****************************************************************/
size_t SubsetVectorData::numberOfDataBlocks() const { return d_dataBlockSize.size(); }
size_t SubsetVectorData::sizeOfDataBlock( size_t i ) const { return d_dataBlockSize[i]; }
void *SubsetVectorData::getRawDataBlockAsVoid( size_t i )
{
    double *ptr = d_dataBlockPtr[i];
    return (void *) ptr;
}
const void *SubsetVectorData::getRawDataBlockAsVoid( size_t i ) const
{
    double *ptr = d_dataBlockPtr[i];
    return (const void *) ptr;
}
void SubsetVectorData::putRawData( const void *in, const typeID &id )
{
    AMP_ASSERT( id == getTypeID<double>() );
    auto data = reinterpret_cast<const double *>( in );
    size_t k  = 0;
    for ( size_t i = 0; i < d_dataBlockPtr.size(); i++ ) {
        for ( size_t j = 0; j < d_dataBlockSize[i]; j++, k++ )
            d_dataBlockPtr[i][j] = data[k];
    }
}
void SubsetVectorData::copyOutRawData( void *out, const typeID &id ) const
{
    AMP_ASSERT( id == getTypeID<double>() );
    auto data = reinterpret_cast<double *>( out );
    size_t k  = 0;
    for ( size_t i = 0; i < d_dataBlockPtr.size(); i++ ) {
        for ( size_t j = 0; j < d_dataBlockSize[i]; j++, k++ )
            data[k] = d_dataBlockPtr[i][j];
    }
}


/****************************************************************
 * Functions get/set/add values                                  *
 ****************************************************************/
void SubsetVectorData::addValuesByLocalID( size_t N, const size_t *ndx, const double *vals )
{
    constexpr size_t N_max = 128;
    while ( N > N_max ) {
        addValuesByLocalID( N_max, ndx, vals );
        N -= N_max;
        ndx  = &ndx[N_max];
        vals = &vals[N_max];
    }
    AMP_ASSERT( d_ViewVector );
    size_t index[N_max];
    for ( size_t i = 0; i != N; i++ )
        index[i] = d_SubsetLocalIDToViewGlobalID[ndx[i]] - d_parentLocalStartID;
    d_ViewVector->addValuesByLocalID( N, index, vals );
}
void SubsetVectorData::setValuesByLocalID( size_t N, const size_t *ndx, const double *vals )
{
    constexpr size_t N_max = 128;
    while ( N > N_max ) {
        setValuesByLocalID( N_max, ndx, vals );
        N -= N_max;
        ndx  = &ndx[N_max];
        vals = &vals[N_max];
    }
    AMP_ASSERT( d_ViewVector );
    size_t index[N_max];
    for ( size_t i = 0; i != N; i++ )
        index[i] = d_SubsetLocalIDToViewGlobalID[ndx[i]] - d_parentLocalStartID;
    d_ViewVector->setValuesByLocalID( N, index, vals );
}
void SubsetVectorData::getValuesByLocalID( size_t N, const size_t *ndx, double *vals ) const
{
    constexpr size_t N_max = 128;
    while ( N > N_max ) {
        getValuesByLocalID( N_max, ndx, vals );
        N -= N_max;
        ndx  = &ndx[N_max];
        vals = &vals[N_max];
    }
    AMP_ASSERT( d_ViewVector );
    size_t index[N_max];
    for ( size_t i = 0; i != N; i++ )
        index[i] = d_SubsetLocalIDToViewGlobalID[ndx[i]] - d_parentLocalStartID;
    d_ViewVector->getValuesByLocalID( N, index, vals );
}

void SubsetVectorData::swapData( VectorData &rhs )
{
    auto s = dynamic_cast<SubsetVectorData *>( &rhs );
    AMP_ASSERT( s != nullptr );
    std::swap( d_ViewVector, s->d_ViewVector );
    std::swap( d_SubsetLocalIDToViewGlobalID, s->d_SubsetLocalIDToViewGlobalID );
    std::swap( d_dataBlockSize, s->d_dataBlockSize );
    std::swap( d_dataBlockPtr, s->d_dataBlockPtr );
}
std::shared_ptr<VectorData> SubsetVectorData::cloneData() const
{
    AMP_ERROR( "Not finished" );
    return std::shared_ptr<VectorData>();
}


} // namespace AMP::LinearAlgebra
