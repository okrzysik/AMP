#include "AMP/operators/map/AsyncMapOperator.h"
#include "AMP/mesh/MultiMesh.h"
#include "AMP/operators/map/AsyncMapOperatorParameters.h"

#include "ProfilerApp.h"


namespace AMP::Operator {


AsyncMapOperator::AsyncMapOperator( std::shared_ptr<const OperatorParameters> p )
    : AsynchronousOperator( p )
{
    // Fill some basic info
    auto params = std::dynamic_pointer_cast<const AsyncMapOperatorParameters>( p );
    d_MapComm   = params->d_MapComm;
    d_mesh1     = params->d_Mesh1;
    d_mesh2     = params->d_Mesh2;
    AMP_INSIST( !d_MapComm.isNull(), "NULL communicator for map is invalid" );
    AMP_INSIST( d_MapComm.anyReduce( bool( d_mesh1 ) ), "Somebody must own mesh 1" );
    AMP_INSIST( d_MapComm.anyReduce( bool( d_mesh2 ) ), "Somebody must own mesh 2" );
    // Create a multimesh to use for the operator base class for subsetting
    std::vector<std::shared_ptr<AMP::Mesh::Mesh>> meshes;
    if ( d_mesh1 )
        meshes.push_back( d_mesh1 );
    if ( d_mesh2 )
        meshes.push_back( d_mesh2 );
    d_Mesh = std::make_shared<AMP::Mesh::MultiMesh>( "mesh", d_MapComm, meshes );
    // Get the input variable
    bool var  = params->d_db->keyExists( "VariableName" );
    bool var1 = params->d_db->keyExists( "VariableName1" );
    bool var2 = params->d_db->keyExists( "VariableName2" );
    AMP_INSIST( var1 || var2 || var, "VariableName must exist in database" );
    if ( var ) {
        AMP_INSIST( !var1 && !var2,
                    "VariableName is used, VariableName1 and VariableName2cannot be used" );
        std::string variableName = params->d_db->getString( "VariableName" );
        d_inpVariable            = std::make_shared<AMP::LinearAlgebra::Variable>( variableName );
        d_outVariable            = std::make_shared<AMP::LinearAlgebra::Variable>( variableName );
    } else {
        AMP_INSIST( var1 && var2, "Both VariableName1 and VariableName2 must be used" );
        std::string variableName1 = params->d_db->getString( "VariableName1" );
        std::string variableName2 = params->d_db->getString( "VariableName2" );
        d_inpVariable             = std::make_shared<AMP::LinearAlgebra::Variable>( variableName1 );
        d_outVariable             = std::make_shared<AMP::LinearAlgebra::Variable>( variableName2 );
    }
}


AsyncMapOperator::~AsyncMapOperator() = default;


void AsyncMapOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                              AMP::LinearAlgebra::Vector::shared_ptr f )
{
    PROFILE_START( "apply" );
    applyStart( u, f );
    applyFinish( u, f );
    if ( requiresMakeConsistentSet() ) {
        AMP_ASSERT( d_OutputVector );
        d_OutputVector->makeConsistent(
            AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    }
    PROFILE_STOP( "apply" );
}


bool AsyncMapOperator::requiresMakeConsistentSet() { return false; }

std::shared_ptr<AMP::Mesh::Mesh> AsyncMapOperator::getMesh( int which )
{
    if ( which == 1 ) {
        return d_mesh1;
    } else if ( which == 2 ) {
        return d_mesh2;
    } else {
        AMP_ERROR( "Wrong option!" );
        return std::shared_ptr<AMP::Mesh::Mesh>();
    }
}
} // namespace AMP::Operator
