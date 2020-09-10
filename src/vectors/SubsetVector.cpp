#include <algorithm>

#include "AMP/discretization/subsetDOFManager.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/SubsetVariable.h"
#include "AMP/vectors/SubsetVector.h"
#include "AMP/vectors/SubsetVectorData.h"
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
    std::shared_ptr<SubsetVector> retVal( new SubsetVector() );
    retVal->setVariable( var );
    auto params          = std::make_shared<SubsetVectorParameters>();
    params->d_ViewVector = std::const_pointer_cast<Vector>( v );
    params->d_DOFManager = subsetDOF;
    params->d_CommList   = commList;
    retVal->d_VectorData = std::make_shared<SubsetVectorData>( params );
    retVal->d_DOFManager = subsetDOF; // at present all vectors need to have a dof manager

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


std::string SubsetVector::type() const
{
    std::string retVal = "Subset Vector";
    retVal += " ( view of ";
    retVal += d_VectorData->VectorDataName();
    retVal += " )";
    return retVal;
}


void SubsetVector::swapVectors( Vector &rhs )
{
#if 1
    d_VectorData->swapData( *( rhs.getVectorData() ) );
#else
    auto s = dynamic_cast<SubsetVector *>( &rhs );
    AMP_ASSERT( s != nullptr );
    std::swap( d_ViewVector, s->d_ViewVector );
    std::swap( d_SubsetLocalIDToViewGlobalID, s->d_SubsetLocalIDToViewGlobalID );
    std::swap( d_dataBlockSize, s->d_dataBlockSize );
    std::swap( d_dataBlockPtr, s->d_dataBlockPtr );
#endif
}


} // namespace LinearAlgebra
} // namespace AMP
