#include "AMP/operators/Operator.h"
#include "AMP/utils/Utilities.h"
#include "ProfilerApp.h"


namespace AMP {
namespace Operator {


int Operator::d_iInstance_id = 0;


Operator::Operator()
{
    d_iObject_id           = Operator::d_iInstance_id;
    d_iDebugPrintInfoLevel = 0;
    Operator::d_iInstance_id++;
}


Operator::Operator( const std::shared_ptr<OperatorParameters> &params )
{
    AMP_INSIST( ( ( params.get() ) != nullptr ), "NULL parameter" );

    d_Mesh                 = params->d_Mesh;
    d_iObject_id           = Operator::d_iInstance_id;
    d_iDebugPrintInfoLevel = 0;
    Operator::d_iInstance_id++;

    // try and keep the next call the last in the function
    // so as not to override any parameters set through it
    // by accident
    getFromInput( params->d_db );
}


void Operator::reset( const std::shared_ptr<OperatorParameters> &params )
{
    AMP_INSIST( ( ( params.get() ) != nullptr ), "NULL parameter" );

    // try and keep the next call the last in the function
    // so as not to override any parameters set through it
    // by accident
    getFromInput( params->d_db );
}


void Operator::residual( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                         AMP::LinearAlgebra::Vector::const_shared_ptr u,
                         AMP::LinearAlgebra::Vector::shared_ptr r )
{
    AMP_INSIST( ( ( u.get() ) != nullptr ), "NULL Solution Vector" );
    AMP_INSIST( ( ( r.get() ) != nullptr ), "NULL Residual Vector" );

    apply( u, r );

    auto rInternal = subsetOutputVector( r );
    AMP_INSIST( ( rInternal.get() != nullptr ), "rInternal is NULL" );

    // the rhs can be NULL
    if ( f.get() != nullptr ) {
        auto fInternal = subsetOutputVector( f );
        rInternal->subtract( fInternal, rInternal, rInternal );
    } else {
        rInternal->scale( -1.0 );
    }

    rInternal->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
}

void Operator::getFromInput( std::shared_ptr<AMP::Database> db )
{
    AMP_INSIST( ( ( db.get() ) != nullptr ), "NULL database" );

    d_iDebugPrintInfoLevel = db->getWithDefault( "print_info_level", 0 );
}


AMP::LinearAlgebra::Vector::shared_ptr
Operator::subsetOutputVector( AMP::LinearAlgebra::Vector::shared_ptr vec )
{
    PROFILE_START( "subsetOutputVector", 1 );
    auto var = getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr retvec;
    if ( d_Mesh.get() != nullptr ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        auto meshSubsetVec = vec->select( meshSelector, vec->getVariable()->getName() );
        auto varSubsetVec  = meshSubsetVec->subsetVectorForVariable( var );
        retvec             = varSubsetVec;
    } else {
        retvec = vec->subsetVectorForVariable( var );
    }
    PROFILE_STOP( "subsetOutputVector", 1 );
    return retvec;
}


AMP::LinearAlgebra::Vector::shared_ptr
Operator::subsetInputVector( AMP::LinearAlgebra::Vector::shared_ptr vec )
{
    PROFILE_START( "subsetInputVector", 1 );
    auto var = getInputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr retvec;
    if ( d_Mesh.get() != nullptr ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        auto meshSubsetVec = vec->select( meshSelector, vec->getVariable()->getName() );
        auto varSubsetVec  = meshSubsetVec->subsetVectorForVariable( var );
        retvec             = varSubsetVec;
    } else {
        retvec = vec->subsetVectorForVariable( var );
    }
    PROFILE_STOP( "subsetInputVector", 1 );
    return retvec;
}


AMP::LinearAlgebra::Vector::const_shared_ptr
Operator::subsetOutputVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec )
{
    PROFILE_START( "constSubsetOutputVector", 1 );
    auto var = getOutputVariable();
    AMP::LinearAlgebra::Vector::const_shared_ptr retvec;
    if ( d_Mesh.get() != nullptr ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        auto meshSubsetVec = vec->constSelect( meshSelector, vec->getVariable()->getName() );
        auto varSubsetVec  = meshSubsetVec->constSubsetVectorForVariable( var );
        retvec             = varSubsetVec;
    } else {
        retvec = vec->constSubsetVectorForVariable( var );
    }
    PROFILE_STOP( "constSubsetOutputVector", 1 );
    return retvec;
}


AMP::LinearAlgebra::Vector::const_shared_ptr
Operator::subsetInputVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec )
{
    PROFILE_START( "constSubsetInputVector", 1 );
    auto var = getInputVariable();
    AMP::LinearAlgebra::Vector::const_shared_ptr retvec;
    if ( d_Mesh.get() != nullptr ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        auto meshSubsetVec = vec->constSelect( meshSelector, vec->getVariable()->getName() );
        auto varSubsetVec  = meshSubsetVec->constSubsetVectorForVariable( var );
        retvec             = varSubsetVec;
    } else {
        retvec = vec->constSubsetVectorForVariable( var );
    }
    PROFILE_STOP( "constSubsetInputVector", 1 );
    return retvec;
}
} // namespace Operator
} // namespace AMP
