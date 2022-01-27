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
    auto parentDOF = d_ViewVector->getDOFManager();
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
void SubsetVectorData::putRawData( const double *in )
{
    size_t k = 0;
    for ( size_t i = 0; i < d_dataBlockPtr.size(); i++ ) {
        for ( size_t j = 0; j < d_dataBlockSize[i]; j++, k++ )
            d_dataBlockPtr[i][j] = in[k];
    }
}
void SubsetVectorData::copyOutRawData( double *out ) const
{
    size_t k = 0;
    for ( size_t i = 0; i < d_dataBlockPtr.size(); i++ ) {
        for ( size_t j = 0; j < d_dataBlockSize[i]; j++, k++ )
            out[k] = d_dataBlockPtr[i][j];
    }
}

/****************************************************************
 * Functions add/set values by ID                                *
 ****************************************************************/
void SubsetVectorData::addLocalValuesByGlobalID( int cnt, size_t *ndx, const double *vals )
{
    AMP_ASSERT( d_ViewVector );
    if ( cnt == 0 )
        return;
    std::vector<size_t> parentDOFs;
    if ( d_SubsetLocalIDToViewGlobalID.size() == d_ViewVector->getLocalSize() ) {
        parentDOFs.resize( cnt );
        int offset = d_ViewVector->getDOFManager()->beginDOF() - d_DOFManager->beginDOF();
        for ( int i = 0; i < cnt; i++ )
            parentDOFs[i] = ndx[i] + offset;
    } else {
        auto DOFManager =
            std::dynamic_pointer_cast<AMP::Discretization::subsetDOFManager>( d_DOFManager );
        std::vector<size_t> subsetDOFs( cnt );
        for ( int i = 0; i < cnt; i++ )
            subsetDOFs[i] = ndx[i];
        parentDOFs = DOFManager->getParentDOF( subsetDOFs );
    }
    d_ViewVector->addLocalValuesByGlobalID( cnt, &parentDOFs[0], vals );
}
void SubsetVectorData::getLocalValuesByGlobalID( int cnt, size_t *ndx, double *vals ) const
{
    AMP_ASSERT( d_ViewVector );
    if ( cnt == 0 )
        return;
    std::vector<size_t> parentDOFs;
    if ( d_SubsetLocalIDToViewGlobalID.size() == d_ViewVector->getLocalSize() ) {
        parentDOFs.resize( cnt );
        int offset = d_ViewVector->getDOFManager()->beginDOF() - d_DOFManager->beginDOF();
        for ( int i = 0; i < cnt; i++ )
            parentDOFs[i] = ndx[i] + offset;
    } else {
        auto DOFManager =
            std::dynamic_pointer_cast<AMP::Discretization::subsetDOFManager>( d_DOFManager );
        std::vector<size_t> subsetDOFs( cnt );
        for ( int i = 0; i < cnt; i++ )
            subsetDOFs[i] = ndx[i];
        parentDOFs = DOFManager->getParentDOF( subsetDOFs );
    }
    d_ViewVector->getLocalValuesByGlobalID( cnt, &parentDOFs[0], vals );
}
void SubsetVectorData::setLocalValuesByGlobalID( int cnt, size_t *ndx, const double *vals )
{
    AMP_ASSERT( d_ViewVector );
    if ( cnt == 0 )
        return;
    std::vector<size_t> parentDOFs;
    if ( d_SubsetLocalIDToViewGlobalID.size() == d_ViewVector->getLocalSize() ) {
        parentDOFs.resize( cnt );
        int offset = d_ViewVector->getDOFManager()->beginDOF() - d_DOFManager->beginDOF();
        for ( int i = 0; i < cnt; i++ )
            parentDOFs[i] = ndx[i] + offset;
    } else {
        auto DOFManager =
            std::dynamic_pointer_cast<AMP::Discretization::subsetDOFManager>( d_DOFManager );
        std::vector<size_t> subsetDOFs( cnt );
        for ( int i = 0; i < cnt; i++ )
            subsetDOFs[i] = ndx[i];
        parentDOFs = DOFManager->getParentDOF( subsetDOFs );
    }
    d_ViewVector->setLocalValuesByGlobalID( cnt, &parentDOFs[0], vals );
}
void SubsetVectorData::addValuesByLocalID( int cnt, size_t *ndx, const double *vals )
{
    AMP_ASSERT( d_ViewVector );
    auto t = new size_t[cnt];
    for ( int i = 0; i != cnt; i++ )
        t[i] = d_SubsetLocalIDToViewGlobalID[ndx[i]];
    d_ViewVector->addValuesByLocalID( cnt, t, vals );
    delete[] t;
}
void SubsetVectorData::setValuesByLocalID( int cnt, size_t *ndx, const double *vals )
{
    AMP_ASSERT( d_ViewVector );
    auto t = new size_t[cnt];
    for ( int i = 0; i != cnt; i++ )
        t[i] = d_SubsetLocalIDToViewGlobalID[ndx[i]];
    d_ViewVector->setValuesByLocalID( cnt, t, vals );
    delete[] t;
}
void SubsetVectorData::getValuesByLocalID( int cnt, size_t *ndx, double *vals ) const
{
    AMP_ASSERT( d_ViewVector );
    auto t = new size_t[cnt];
    for ( int i = 0; i != cnt; i++ )
        t[i] = d_SubsetLocalIDToViewGlobalID[ndx[i]];
    d_ViewVector->getValuesByLocalID( cnt, t, vals );
    delete[] t;
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
