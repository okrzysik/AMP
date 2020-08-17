#include <algorithm>

#include "AMP/discretization/subsetDOFManager.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/SubsetVariable.h"
#include "AMP/vectors/SubsetVector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/data/VectorDataIterator.h"

#include "ProfilerApp.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Contructors                                                   *
 ****************************************************************/
Vector::shared_ptr SubsetVector::view( Vector::shared_ptr v, Variable::shared_ptr var_in )
{
    return std::const_pointer_cast<Vector>(
        SubsetVector::view( Vector::const_shared_ptr( v ), var_in ) );
}
Vector::const_shared_ptr SubsetVector::view( Vector::const_shared_ptr v,
                                             Variable::shared_ptr var_in )
{
    PROFILE_START( "view", 2 );
    auto var = std::dynamic_pointer_cast<SubsetVariable>( var_in );
    AMP_ASSERT( var.get() != nullptr );
    if ( std::dynamic_pointer_cast<const MultiVector>( v ) ) {
        // We are dealing with a multivector, it is more efficient to subset the individual pieces,
        // then create a new mulitvector
        auto mv = std::dynamic_pointer_cast<const MultiVector>( v );
        std::vector<Vector::const_shared_ptr> vec_list;
        for ( size_t i = 0; i < mv->getNumberOfSubvectors(); i++ ) {
            Vector::const_shared_ptr sub_vec = SubsetVector::view( mv->getVector( i ), var_in );
            if ( sub_vec != nullptr )
                vec_list.push_back( sub_vec );
        }
        if ( vec_list.empty() ) {
            PROFILE_STOP2( "view", 2 );
            return Vector::const_shared_ptr();
        }
        auto parentDOF = v->getDOFManager();
        auto subsetDOF = var->getSubsetDOF( parentDOF );
        PROFILE_STOP2( "view", 2 );
        return MultiVector::const_create(
            v->getVariable()->getName(), subsetDOF->getComm(), vec_list );
    }
    // Subset the DOFManager and create a new communication list
    auto parentDOF = v->getDOFManager();
    auto subsetDOF = var->getSubsetDOF( parentDOF );
    if ( subsetDOF.get() == nullptr ) {
        PROFILE_STOP2( "view", 2 );
        return Vector::shared_ptr();
    } else if ( subsetDOF->numGlobalDOF() == 0 ) {
        PROFILE_STOP2( "view", 2 );
        return Vector::shared_ptr();
    } else if ( subsetDOF == parentDOF ) {
        PROFILE_STOP2( "view", 2 );
        return v;
    }
    auto remote_DOFs = subsetDOF->getRemoteDOFs();
    bool ghosts      = subsetDOF->getComm().anyReduce( !remote_DOFs.empty() );
    AMP::LinearAlgebra::CommunicationList::shared_ptr commList;
    if ( !ghosts ) {
        commList = AMP::LinearAlgebra::CommunicationList::createEmpty( subsetDOF->numLocalDOF(),
                                                                       subsetDOF->getComm() );
    } else {
        // Construct the communication list
        auto params           = std::make_shared<AMP::LinearAlgebra::CommunicationListParameters>();
        params->d_comm        = subsetDOF->getComm();
        params->d_localsize   = subsetDOF->numLocalDOF();
        params->d_remote_DOFs = remote_DOFs;
        commList              = std::make_shared<AMP::LinearAlgebra::CommunicationList>( params );
    }
    // Create the new subset vector
    auto retVal = std::make_shared<SubsetVector>();
    retVal->setVariable( var );
    retVal->d_ViewVector = std::const_pointer_cast<Vector>( v );
    retVal->d_DOFManager = subsetDOF;
    retVal->setCommunicationList( commList );
    auto tmp = std::dynamic_pointer_cast<AMP::Discretization::subsetDOFManager>( subsetDOF );
    if ( tmp != nullptr ) {
        retVal->d_SubsetLocalIDToViewGlobalID = tmp->getLocalParentDOFs();
    } else if ( subsetDOF->numLocalDOF() == parentDOF->numLocalDOF() ) {
        retVal->d_SubsetLocalIDToViewGlobalID.resize( parentDOF->numLocalDOF() );
        size_t beginDOF = parentDOF->beginDOF();
        size_t endDOF   = parentDOF->endDOF();
        for ( size_t i = beginDOF; i < endDOF; i++ )
            retVal->d_SubsetLocalIDToViewGlobalID[i - beginDOF] = i;
    } else {
        AMP_ERROR( "Internal error with SubsetVector" );
    }
    // Get a pointer to every value in the subset
    std::vector<double *> data_ptr( retVal->d_SubsetLocalIDToViewGlobalID.size(), nullptr );
    auto iterator   = retVal->d_ViewVector->constBegin();
    size_t last_pos = retVal->d_ViewVector->getCommunicationList()->getStartGID();
    for ( size_t i = 0; i < data_ptr.size(); i++ ) {
        iterator += (int) ( retVal->d_SubsetLocalIDToViewGlobalID[i] - last_pos );
        last_pos    = retVal->d_SubsetLocalIDToViewGlobalID[i];
        data_ptr[i] = const_cast<double *>( &( *iterator ) );
    }
    // Create the data blocks
    // For now use one datablock for each value, this needs to be changed
    retVal->d_dataBlockPtr  = data_ptr;
    retVal->d_dataBlockSize = std::vector<size_t>( data_ptr.size(), 1 );
    // We should decide on a better way to set the vector operations
    // for efficiency
    retVal->d_VectorOps = std::make_shared<VectorOperationsDefault<double>>();
    PROFILE_STOP( "view", 2 );
    return retVal;
}
Vector::shared_ptr SubsetVector::cloneVector( Variable::shared_ptr var ) const
{
    // Ideally this function should create a new dense vector of the same type as d_ViewVector
    // For now, create a dense vector of a possibly new type
    Vector::shared_ptr vec = AMP::LinearAlgebra::createVector( d_DOFManager, var );
    return vec;
}


/****************************************************************
 * Functions to access the raw data blocks                       *
 ****************************************************************/
size_t SubsetVector::numberOfDataBlocks() const { return d_dataBlockSize.size(); }
size_t SubsetVector::sizeOfDataBlock( size_t i ) const { return d_dataBlockSize[i]; }
void *SubsetVector::getRawDataBlockAsVoid( size_t i )
{
    double *ptr = d_dataBlockPtr[i];
    return (void *) ptr;
}
const void *SubsetVector::getRawDataBlockAsVoid( size_t i ) const
{
    double *ptr = d_dataBlockPtr[i];
    return (const void *) ptr;
}
void SubsetVector::putRawData( const double *in )
{
    size_t k = 0;
    for ( size_t i = 0; i < d_dataBlockPtr.size(); i++ ) {
        for ( size_t j = 0; j < d_dataBlockSize[i]; j++, k++ )
            d_dataBlockPtr[i][j] = in[k];
    }
}
void SubsetVector::copyOutRawData( double *out ) const
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
void SubsetVector::addLocalValuesByGlobalID( int cnt, size_t *ndx, const double *vals )
{
    AMP_ASSERT( d_ViewVector.get() != nullptr );
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
void SubsetVector::getLocalValuesByGlobalID( int cnt, size_t *ndx, double *vals ) const
{
    AMP_ASSERT( d_ViewVector.get() != nullptr );
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
void SubsetVector::setLocalValuesByGlobalID( int cnt, size_t *ndx, const double *vals )
{
    AMP_ASSERT( d_ViewVector.get() != nullptr );
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
void SubsetVector::addValuesByLocalID( int cnt, size_t *ndx, const double *vals )
{
    AMP_ASSERT( d_ViewVector.get() != nullptr );
    auto t = new size_t[cnt];
    for ( int i = 0; i != cnt; i++ )
        t[i] = d_SubsetLocalIDToViewGlobalID[ndx[i]];
    d_ViewVector->addValuesByLocalID( cnt, t, vals );
    delete[] t;
}
void SubsetVector::setValuesByLocalID( int cnt, size_t *ndx, const double *vals )
{
    AMP_ASSERT( d_ViewVector.get() != nullptr );
    auto t = new size_t[cnt];
    for ( int i = 0; i != cnt; i++ )
        t[i] = d_SubsetLocalIDToViewGlobalID[ndx[i]];
    d_ViewVector->setValuesByLocalID( cnt, t, vals );
    delete[] t;
}
void SubsetVector::getValuesByLocalID( int cnt, size_t *ndx, double *vals ) const
{
    AMP_ASSERT( d_ViewVector.get() != nullptr );
    auto t = new size_t[cnt];
    for ( int i = 0; i != cnt; i++ )
        t[i] = d_SubsetLocalIDToViewGlobalID[ndx[i]];
    d_ViewVector->getValuesByLocalID( cnt, t, vals );
    delete[] t;
}

size_t SubsetVector::getLocalSize() const { return getCommunicationList()->numLocalRows(); }
size_t SubsetVector::getGlobalSize() const { return getCommunicationList()->getTotalSize(); }


void SubsetVector::aliasVector( Vector & ) { AMP_ERROR( "cannot alias a subset vector.....yet" ); }


std::string SubsetVector::type() const
{
    std::string retVal = "Subset Vector";
    if ( d_ViewVector ) {
        retVal += " ( view of ";
        retVal += d_ViewVector->type();
        retVal += " )";
    }
    return retVal;
}


void SubsetVector::swapVectors( Vector &rhs )
{
    auto s = dynamic_cast<SubsetVector *>( &rhs );
    AMP_ASSERT( s != nullptr );
    std::swap( d_ViewVector, s->d_ViewVector );
    std::swap( d_SubsetLocalIDToViewGlobalID, s->d_SubsetLocalIDToViewGlobalID );
    std::swap( d_dataBlockSize, s->d_dataBlockSize );
    std::swap( d_dataBlockPtr, s->d_dataBlockPtr );
}


} // namespace LinearAlgebra
} // namespace AMP
