#include "AMP/operators/Operator.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorSelector.h"

#include "ProfilerApp.h"


namespace AMP::Operator {


int Operator::d_iInstance_id = 0;


Operator::Operator()
{
    d_iObject_id           = Operator::d_iInstance_id;
    d_iDebugPrintInfoLevel = 0;
    Operator::d_iInstance_id++;
    d_memory_location = AMP::Utilities::MemoryType::none;
}


Operator::Operator( std::shared_ptr<const OperatorParameters> params )
{
    AMP_INSIST( params, "NULL parameter" );

    d_Mesh                 = params->d_Mesh;
    d_iObject_id           = Operator::d_iInstance_id;
    d_iDebugPrintInfoLevel = 0;
    Operator::d_iInstance_id++;

    // try and keep the next call the last in the function
    // so as not to override any parameters set through it
    // by accident
    getFromInput( params->d_db );
}


void Operator::reset( std::shared_ptr<const OperatorParameters> params )
{
    AMP_INSIST( params, "NULL parameter" );

    // try and keep the next call the last in the function
    // so as not to override any parameters set through it
    // by accident
    getFromInput( params->d_db );
}


void Operator::residual( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                         AMP::LinearAlgebra::Vector::const_shared_ptr u,
                         AMP::LinearAlgebra::Vector::shared_ptr r )
{
    AMP_INSIST( u, "NULL Solution Vector" );
    AMP_INSIST( r, "NULL Residual Vector" );
    AMP_ASSERT( u->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::UNCHANGED );

    apply( u, r );

    auto rInternal = subsetOutputVector( r );
    AMP_INSIST( ( rInternal ), "rInternal is NULL" );

    // the rhs can be NULL
    if ( f ) {
        auto fInternal = subsetOutputVector( f );
        rInternal->subtract( *fInternal, *rInternal );
    } else {
        rInternal->scale( -1.0 );
    }

    rInternal->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}


std::shared_ptr<OperatorParameters>
Operator::getParameters( const std::string &type,
                         std::shared_ptr<const AMP::LinearAlgebra::Vector> u,
                         std::shared_ptr<OperatorParameters> )
{
    if ( type == "Jacobian" ) {
        return getJacobianParameters( u );
    } else {
        // Derived class should implement this
        AMP_ERROR( "Unknown OperatorParameters type specified" );
    }
    return nullptr;
}


void Operator::getFromInput( std::shared_ptr<AMP::Database> db )
{
    AMP_INSIST( db, "NULL database" );

    d_iDebugPrintInfoLevel = db->getWithDefault<int>( "print_info_level", 0 );

    if ( d_memory_location == AMP::Utilities::MemoryType::none ) {
        auto memLoc       = db->getWithDefault<std::string>( "MemoryLocation", "host" );
        d_memory_location = AMP::Utilities::memoryLocationFromString( memLoc );
    }
    if ( d_backend == AMP::Utilities::Backend::none ) {
        if ( db->keyExists( "AccelerationBackend" ) ) {
            auto bcknd = db->getString( "AccelerationBackend" );
            d_backend  = AMP::Utilities::backendFromString( bcknd );
        } else {
            d_backend = AMP::Utilities::getDefaultBackend( d_memory_location );
        }
    }
}


AMP::LinearAlgebra::Vector::shared_ptr
Operator::subsetOutputVector( AMP::LinearAlgebra::Vector::shared_ptr vec )
{
    PROFILE( "subsetOutputVector", 1 );
    // Subset for mesh (if set)
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        vec = vec->select( meshSelector );
    }
    // Subset for variable (if set)
    auto var = getOutputVariable();
    if ( var )
        vec = vec->subsetVectorForVariable( var );
    return vec;
}


AMP::LinearAlgebra::Vector::shared_ptr
Operator::subsetInputVector( AMP::LinearAlgebra::Vector::shared_ptr vec )
{
    PROFILE( "subsetInputVector", 1 );
    // Subset for mesh (if set)
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        vec = vec->select( meshSelector );
    }
    // Subset for variable (if set)
    auto var = getInputVariable();
    if ( var )
        vec = vec->subsetVectorForVariable( var );
    return vec;
}


AMP::LinearAlgebra::Vector::const_shared_ptr
Operator::subsetOutputVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec )
{
    PROFILE( "constSubsetOutputVector", 1 );
    // Subset for mesh (if set)
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        vec = vec->select( meshSelector );
    }
    // Subset for variable (if set)
    auto var = getOutputVariable();
    if ( var )
        vec = vec->subsetVectorForVariable( var );
    return vec;
}


AMP::LinearAlgebra::Vector::const_shared_ptr
Operator::subsetInputVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec )
{
    PROFILE( "constSubsetInputVector", 1 );
    // Subset for mesh (if set)
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        vec = vec->select( meshSelector );
    }
    // Subset for variable (if set)
    auto var = getInputVariable();
    if ( var )
        vec = vec->subsetVectorForVariable( var );
    return vec;
}

void Operator::makeConsistent( std::shared_ptr<AMP::LinearAlgebra::Vector> vec )
{
    AMP_ASSERT( vec );
    vec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}

} // namespace AMP::Operator
