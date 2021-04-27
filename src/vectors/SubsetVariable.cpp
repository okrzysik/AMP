#include "AMP/vectors/SubsetVariable.h"
#include "AMP/discretization/subsetDOFManager.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/data/SubsetVectorData.h"
#include "AMP/vectors/data/VectorDataIterator.h"

#include "ProfilerApp.h"

#include <algorithm>


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Contructors                                                   *
 ****************************************************************/
Vector::shared_ptr SubsetVariable::view( Vector::shared_ptr v, Variable::shared_ptr var_in )
{
    return std::const_pointer_cast<Vector>(
        SubsetVariable::view( Vector::const_shared_ptr( v ), var_in ) );
}
Vector::const_shared_ptr SubsetVariable::view( Vector::const_shared_ptr v,
                                               Variable::shared_ptr var_in )
{
    PROFILE_START( "view", 2 );
    auto var = std::dynamic_pointer_cast<SubsetVariable>( var_in );
    AMP_ASSERT( var );
    if ( std::dynamic_pointer_cast<const MultiVector>( v ) ) {
        // We are dealing with a multivector, it is more efficient to subset the individual pieces,
        // then create a new mulitvector
        auto mv = std::dynamic_pointer_cast<const MultiVector>( v );
        std::vector<Vector::const_shared_ptr> vec_list;
        for ( size_t i = 0; i < mv->getNumberOfSubvectors(); i++ ) {
            auto sub_vec = SubsetVariable::view( mv->getVector( i ), var_in );
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
    auto ops             = std::make_shared<VectorOperationsDefault<double>>();
    auto params          = std::make_shared<SubsetVectorParameters>();
    params->d_ViewVector = std::const_pointer_cast<Vector>( v );
    params->d_DOFManager = subsetDOF;
    params->d_CommList   = commList;
    auto data            = std::make_shared<SubsetVectorData>( params );
    auto retVal          = std::make_shared<Vector>( data, ops, var, subsetDOF );
    PROFILE_STOP( "view", 2 );
    return retVal;
}


} // namespace LinearAlgebra
} // namespace AMP
