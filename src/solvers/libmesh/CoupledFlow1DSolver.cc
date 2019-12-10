#include "AMP/solvers/libmesh/CoupledFlow1DSolver.h"
#include "AMP/operators/map/Map1Dto3D.h"
#include "AMP/operators/map/Map3Dto1D.h"
#include "AMP/operators/subchannel/FlowFrapconJacobian.h"
#include "AMP/solvers/libmesh/Flow1DSolver.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVector.h"


namespace AMP {
namespace Solver {

CoupledFlow1DSolver::CoupledFlow1DSolver( AMP::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters ), d_numpoints( 0 )
{
    AMP_ASSERT( parameters.get() != nullptr );

    AMP::shared_ptr<CoupledFlow1DSolverParameters> params =
        AMP::dynamic_pointer_cast<CoupledFlow1DSolverParameters>( parameters );

    std::string flowOutVar = ( ( params->d_pOperator )->getOutputVariable() )->getName();

    d_pOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::CoupledFlowFrapconOperator>( params->d_pOperator );

    d_flow1DSolver         = AMP::dynamic_pointer_cast<Flow1DSolver>( params->d_flow1DSolver );
    std::string flowInpVar = ( d_flow1DSolver->getInputVariable() )->getName();

    AMP::shared_ptr<AMP::Database> tmp_db1( new AMP::Database( "Dummy" ) );
    tmp_db1->putScalar( "BoundaryId", 4 );
    tmp_db1->putScalar( "InputVariable", flowOutVar );
    tmp_db1->putScalar( "OutputVariable", flowInpVar );
    AMP::shared_ptr<AMP::Operator::MapOperatorParameters> mapflowInternal3to1Params(
        new AMP::Operator::MapOperatorParameters( tmp_db1 ) );
    mapflowInternal3to1Params->d_MapMesh = ( d_flow1DSolver->getOperator() )->getMesh();
    mapflowInternal3to1Params->d_MapComm = mapflowInternal3to1Params->d_MapMesh->getComm();
    d_flowInternal3to1.reset( new AMP::Operator::Map3Dto1D( mapflowInternal3to1Params ) );

    AMP::shared_ptr<AMP::Database> tmp_db2( new AMP::Database( "Dummy" ) );
    tmp_db2->putScalar( "BoundaryId", 4 );
    tmp_db2->putScalar( "InputVariable", flowInpVar );
    tmp_db2->putScalar( "OutputVariable", flowOutVar );
    AMP::shared_ptr<AMP::Operator::MapOperatorParameters> mapflowInternal1to3Params(
        new AMP::Operator::MapOperatorParameters( tmp_db2 ) );
    mapflowInternal1to3Params->d_MapMesh = ( d_flow1DSolver->getOperator() )->getMesh();
    mapflowInternal1to3Params->d_MapComm = mapflowInternal1to3Params->d_MapMesh->getComm();
    d_flowInternal1to3.reset( new AMP::Operator::Map1Dto3D( mapflowInternal1to3Params ) );

    ( AMP::dynamic_pointer_cast<AMP::Operator::Map3Dto1D>( d_flowInternal3to1 ) )
        ->setZLocations(
            ( AMP::dynamic_pointer_cast<AMP::Operator::Map1Dto3D>( d_flowInternal1to3 ) )
                ->getZLocations() );

    int d_numpoints = ( AMP::dynamic_pointer_cast<AMP::Operator::Map1Dto3D>( d_flowInternal1to3 ) )
                          ->getNumZlocations();
    d_SimpleVariable.reset( new AMP::LinearAlgebra::Variable( flowInpVar ) );

    d_flowInput = AMP::LinearAlgebra::SimpleVector<double>::create( d_numpoints, d_SimpleVariable );
    d_flowOutput =
        AMP::LinearAlgebra::SimpleVector<double>::create( d_numpoints, d_SimpleVariable );
}

CoupledFlow1DSolver::~CoupledFlow1DSolver() = default;

void CoupledFlow1DSolver::setInitialGuess( AMP::shared_ptr<AMP::LinearAlgebra::Vector> ) {}

void CoupledFlow1DSolver::reset( AMP::shared_ptr<SolverStrategyParameters> )
{

    if ( d_pOperator.get() != nullptr ) {}
}

void CoupledFlow1DSolver::resetOperator(
    const AMP::shared_ptr<AMP::Operator::OperatorParameters> params )
{
    if ( d_pOperator.get() != nullptr ) {
        d_pOperator->reset( params );
    }
}

void CoupledFlow1DSolver::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                                 AMP::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    d_inpVariable = d_flow1DSolver->getOperator()->getInputVariable();
    d_outVariable = d_flowInternal1to3->getOutputVariable();

    d_Sol = u->subsetVectorForVariable( d_outVariable );
    d_Rhs = f->constSubsetVectorForVariable( d_outVariable );
    AMP_ASSERT( d_Rhs->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED );

    ( AMP::dynamic_pointer_cast<AMP::Operator::Map3Dto1D>( d_flowInternal3to1 ) )
        ->setVector( d_flowInput );
    ( AMP::dynamic_pointer_cast<AMP::Operator::Map1Dto3D>( d_flowInternal1to3 ) )
        ->setVector( d_Sol );

    d_flowInternal3to1->apply( d_Rhs, nullVec );
    d_flow1DSolver->solve( d_flowInput, d_flowOutput );
    d_flowInternal1to3->apply( d_flowOutput, nullVec );
}
} // namespace Solver
} // namespace AMP
