#include "AMP/operators/subchannel/CoupledFlowFrapconOperator.h"
#include "AMP/operators/subchannel/CoupledFlowFrapconOperatorParameters.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/SimpleVector.h"


namespace AMP {
namespace Operator {

CoupledFlowFrapconOperator::CoupledFlowFrapconOperator(
    const AMP::shared_ptr<OperatorParameters> &params )
    : ColumnOperator( params )
{
    d_Mesh = params->d_Mesh;
    AMP::shared_ptr<CoupledFlowFrapconOperatorParameters> myparams =
        AMP::dynamic_pointer_cast<CoupledFlowFrapconOperatorParameters>( params );
    d_Operators.push_back( myparams->d_Map3to1 );

    std::string flowOutVar = ( ( myparams->d_Map1to3 )->getOutputVariable() )->getName();

    std::string flowInpVar = ( ( myparams->d_FlowOperator )->getOutputVariable() )->getName();
    d_SimpleVariable.reset( new AMP::LinearAlgebra::Variable( flowInpVar ) );

    d_numpoints = ( AMP::dynamic_pointer_cast<AMP::Operator::Map1Dto3D>( myparams->d_Map1to3 ) )
                      ->getNumZlocations();
    d_zPoints.resize( d_numpoints );

    d_flowInput = AMP::LinearAlgebra::SimpleVector<double>::create( d_numpoints, d_SimpleVariable );
    d_flowOutput =
        AMP::LinearAlgebra::SimpleVector<double>::create( d_numpoints, d_SimpleVariable );

    AMP::shared_ptr<AMP::Database> tmp_db1( new AMP::Database( "Dummy" ) );
    tmp_db1->putScalar( "BoundaryId", 4 );
    tmp_db1->putScalar( "InputVariable", flowOutVar );
    tmp_db1->putScalar( "OutputVariable", "FlowInternal" );
    AMP::shared_ptr<AMP::Operator::MapOperatorParameters> mapflowInternal3to1Params(
        new AMP::Operator::MapOperatorParameters( tmp_db1 ) );
    mapflowInternal3to1Params->d_Mesh =
        ( AMP::dynamic_pointer_cast<AMP::Operator::Map3Dto1D>( myparams->d_Map3to1 ) )->getMesh();
    mapflowInternal3to1Params->d_MapMesh =
        ( AMP::dynamic_pointer_cast<AMP::Operator::Map3Dto1D>( myparams->d_Map3to1 ) )->getMesh();
    mapflowInternal3to1Params->d_MapComm = mapflowInternal3to1Params->d_MapMesh->getComm();
    d_flowInternal3to1.reset( new AMP::Operator::Map3Dto1D( mapflowInternal3to1Params ) );

    ( AMP::dynamic_pointer_cast<AMP::Operator::Map3Dto1D>( d_flowInternal3to1 ) )
        ->setVector( d_flowInput );

    d_Operators.push_back( d_flowInternal3to1 );
    d_Operators.push_back( myparams->d_FlowOperator );

    AMP::shared_ptr<AMP::Database> tmp_db2( new AMP::Database( "Dummy" ) );
    tmp_db2->putScalar( "BoundaryId", 4 );
    tmp_db2->putScalar( "InputVariable", "FlowInternal" );
    tmp_db2->putScalar( "OutputVariable", flowOutVar );
    AMP::shared_ptr<AMP::Operator::MapOperatorParameters> mapflowInternal1to3Params(
        new AMP::Operator::MapOperatorParameters( tmp_db2 ) );
    mapflowInternal1to3Params->d_Mesh =
        ( AMP::dynamic_pointer_cast<AMP::Operator::Map3Dto1D>( myparams->d_Map3to1 ) )->getMesh();
    mapflowInternal1to3Params->d_MapMesh =
        ( AMP::dynamic_pointer_cast<AMP::Operator::Map3Dto1D>( myparams->d_Map3to1 ) )->getMesh();
    mapflowInternal1to3Params->d_MapComm = mapflowInternal1to3Params->d_MapMesh->getComm();
    d_flowInternal1to3.reset( new AMP::Operator::Map1Dto3D( mapflowInternal1to3Params ) );

    ( AMP::dynamic_pointer_cast<AMP::Operator::Map3Dto1D>( d_flowInternal3to1 ) )
        ->setZLocations(
            ( AMP::dynamic_pointer_cast<AMP::Operator::Map1Dto3D>( d_flowInternal1to3 ) )
                ->getZLocations() );

    d_Operators.push_back( d_flowInternal1to3 );
    d_Operators.push_back( myparams->d_Map1to3 );
}

void CoupledFlowFrapconOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                        AMP::LinearAlgebra::Vector::shared_ptr r )
{
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    // AMP::LinearAlgebra::Variable::shared_ptr inpVar =
    // (AMP::dynamic_pointer_cast<AMP::Operator::Map1Dto3D>
    // (d_Operators[3]) )->getOutputVariable();
    // AMP::LinearAlgebra::Vector::const_shared_ptr uInternal = subsetInputVector( u );
    AMP::LinearAlgebra::Vector::shared_ptr rInternal = subsetInputVector( r );
    //(AMP::dynamic_pointer_cast<AMP::Operator::Map1Dto3D> (d_Operators[3]))->setVector(uInternal);
    //// Is this
    // necessary
    ( AMP::dynamic_pointer_cast<AMP::Operator::Map1Dto3D>( d_Operators[4] ) )
        ->setVector( rInternal );

    d_Operators[0]->apply( u, nullVec );
    d_Operators[1]->apply( u, nullVec );
    d_Operators[2]->apply( d_flowInput, d_flowOutput );
    // d_Operators[3]->apply(nullVec, d_flowInput, nullVec, a, b);  // Is this necessary
    d_Operators[4]->apply( d_flowOutput, nullVec );
}
} // namespace Operator
} // namespace AMP
