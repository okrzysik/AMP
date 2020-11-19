#include "AMP/solvers/libmesh/PelletStackHelpers.h"

#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/operators/CoupledOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/libmesh/PelletStackOperator.h"
#include "AMP/operators/map/AsyncMapColumnOperator.h"
#include "AMP/operators/map/NodeToNodeMap.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "AMP/vectors/VectorBuilder.h"

#ifdef USE_EXT_PETSC
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#endif

#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"


void helperCreateStackOperatorForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    std::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps,
    std::shared_ptr<AMP::Database> global_input_db,
    std::shared_ptr<AMP::Operator::PelletStackOperator> &pelletStackOp )
{
    auto pelletStackOp_db        = global_input_db->getDatabase( "PelletStackOperator" );
    auto pelletMeshes            = manager->Subset( "PelletMeshes" );
    auto pelletMeshIDs           = pelletMeshes->getBaseMeshIDs();
    unsigned int totalNumPellets = pelletMeshIDs.size();
    pelletStackOp_db->putScalar( "TOTAL_NUMBER_OF_PELLETS", totalNumPellets );
    auto pelletStackOpParams =
        std::make_shared<AMP::Operator::PelletStackOperatorParameters>( pelletStackOp_db );
    pelletStackOpParams->d_pelletStackComm = pelletMeshes->getComm();
    pelletStackOpParams->d_n2nMaps         = n2nmaps;
    pelletStackOpParams->d_Mesh            = pelletMeshes;
    pelletStackOp.reset( new AMP::Operator::PelletStackOperator( pelletStackOpParams ) );
}

void helperCreateColumnOperatorsForPelletMechanics(
    std::vector<unsigned int> localPelletIds,
    std::vector<AMP::Mesh::Mesh::shared_ptr> localMeshes,
    std::shared_ptr<AMP::Database> global_input_db,
    std::shared_ptr<AMP::Operator::ColumnOperator> &nonlinearColumnOperator,
    std::shared_ptr<AMP::Operator::ColumnOperator> &linearColumnOperator )
{
    std::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
    nonlinearColumnOperator.reset( new AMP::Operator::ColumnOperator( emptyParams ) );
    linearColumnOperator.reset( new AMP::Operator::ColumnOperator( emptyParams ) );
    for ( unsigned int id = 0; id < localPelletIds.size(); id++ ) {
        std::string prefix = "";
        if ( localPelletIds[id] == 0 ) {
            prefix = "Bottom";
        }

        auto meshAdapter = localMeshes[id];

        std::shared_ptr<AMP::Operator::ElementPhysicsModel> mechModel;
        auto nonlinOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter,
                prefix + "PelletMechanicsNonlinearBVPOperator",
                global_input_db,
                mechModel ) );
        nonlinearColumnOperator->append( nonlinOperator );

        auto linOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator( meshAdapter,
                                                            prefix +
                                                                "PelletMechanicsLinearBVPOperator",
                                                            global_input_db,
                                                            mechModel ) );
        linearColumnOperator->append( linOperator );
    } // end for id
}

void helperCreateCoupledOperatorForPelletMechanics(
    std::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps,
    std::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator,
    std::shared_ptr<AMP::Operator::CoupledOperator> &coupledOp )
{
    std::shared_ptr<AMP::Database> emptyDb;
    auto coupledOpParams = std::make_shared<AMP::Operator::CoupledOperatorParameters>( emptyDb );
    coupledOpParams->d_MapOperator = n2nmaps;
    coupledOpParams->d_BVPOperator = nonlinearColumnOperator;
    coupledOp.reset( new AMP::Operator::CoupledOperator( coupledOpParams ) );
}

void helperSetFrozenVectorForMapsForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager, std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp )
{
    auto n2nmaps = std::dynamic_pointer_cast<AMP::Operator::AsyncMapColumnOperator>(
        coupledOp->getOperator( 2 ) );
    auto nonlinearColumnOperator =
        std::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    auto dispVar         = nonlinearColumnOperator->getOutputVariable();
    auto nodal3VectorDOF = AMP::Discretization::simpleDOFManager::create(
        manager, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    auto dirichletValues = AMP::LinearAlgebra::createVector( nodal3VectorDOF, dispVar, true );
    if ( n2nmaps ) {
        n2nmaps->setVector( dirichletValues );
    }
    for ( unsigned int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++ ) {
        auto dirichletOp = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
                nonlinearColumnOperator->getOperator( id ) )
                ->getBoundaryOperator() );
        dirichletOp->setDirichletValues( dirichletValues );
    } // end for id
}

void helperCreateAllOperatorsForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    AMP::AMP_MPI,
    std::shared_ptr<AMP::Database> global_input_db,
    std::shared_ptr<AMP::Operator::CoupledOperator> &coupledOp,
    std::shared_ptr<AMP::Operator::ColumnOperator> &linearColumnOperator,
    std::shared_ptr<AMP::Operator::PelletStackOperator> &pelletStackOp )
{
    std::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps;
    if ( global_input_db->keyExists( "MechanicsNodeToNodeMaps" ) ) {
        auto map_db = global_input_db->getDatabase( "MechanicsNodeToNodeMaps" );
        n2nmaps     = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap>(
            manager, map_db );
    }
    helperCreateStackOperatorForPelletMechanics( manager, n2nmaps, global_input_db, pelletStackOp );
    auto localPelletIds = pelletStackOp->getLocalPelletIds();
    auto localMeshes    = pelletStackOp->getLocalMeshes();
    std::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator;
    helperCreateColumnOperatorsForPelletMechanics( localPelletIds,
                                                   localMeshes,
                                                   global_input_db,
                                                   nonlinearColumnOperator,
                                                   linearColumnOperator );
    helperCreateCoupledOperatorForPelletMechanics( n2nmaps, nonlinearColumnOperator, coupledOp );
    helperSetFrozenVectorForMapsForPelletMechanics( manager, coupledOp );
}

void helperCreateVectorsForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr &solVec,
    AMP::LinearAlgebra::Vector::shared_ptr &rhsVec,
    AMP::LinearAlgebra::Vector::shared_ptr &scaledRhsVec )
{
    auto nonlinearColumnOperator =
        std::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    auto dispVar         = nonlinearColumnOperator->getOutputVariable();
    auto nodal3VectorDOF = AMP::Discretization::simpleDOFManager::create(
        manager, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    solVec       = AMP::LinearAlgebra::createVector( nodal3VectorDOF, dispVar, true );
    rhsVec       = AMP::LinearAlgebra::createVector( nodal3VectorDOF, dispVar, true );
    scaledRhsVec = rhsVec->cloneVector();
}

void helperBuildPointLoadRHSForPelletMechanics(
    std::shared_ptr<AMP::Database> global_input_db,
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec )
{
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto nonlinearColumnOperator =
        std::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    rhsVec->zero();
    for ( unsigned int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++ ) {
        auto currOp      = nonlinearColumnOperator->getOperator( id );
        auto meshAdapter = currOp->getMesh();
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
        auto loadOp = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "PointLoad", global_input_db, dummyModel ) );
        loadOp->setVariable( currOp->getOutputVariable() );
        loadOp->apply( nullVec, rhsVec );
    } // end for id
}

void helperApplyBoundaryCorrectionsForPelletMechanics(
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr solVec,
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec )
{
    auto nonlinearColumnOperator =
        std::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    for ( unsigned int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++ ) {
        auto nonlinOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            nonlinearColumnOperator->getOperator( id ) );
        nonlinOperator->modifyInitialSolutionVector( solVec );
        nonlinOperator->modifyRHSvector( rhsVec );
    } // end for id
}

void helperCreateTemperatureVectorsForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    AMP::LinearAlgebra::Vector::shared_ptr &initialTemperatureVec,
    AMP::LinearAlgebra::Vector::shared_ptr &finalTemperatureVec )
{
    auto tempVar        = std::make_shared<AMP::LinearAlgebra::Variable>( "temperature" );
    auto nodalScalarDOF = AMP::Discretization::simpleDOFManager::create(
        manager, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    initialTemperatureVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, tempVar, true );
    finalTemperatureVec   = initialTemperatureVec->cloneVector();
}

void helperSetReferenceTemperatureForPelletMechanics(
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr initialTemperatureVec )
{
    auto nonlinearColumnOperator =
        std::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    for ( unsigned int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++ ) {
        auto bvpOp = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            nonlinearColumnOperator->getOperator( id ) );
        auto mechOp = std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            bvpOp->getVolumeOperator() );
        mechOp->setReferenceTemperature( initialTemperatureVec );
    } // end for id
}

void helperSetFinalTemperatureForPelletMechanics(
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr finalTemperatureVec )
{
    auto nonlinearColumnOperator =
        std::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    for ( unsigned int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++ ) {
        auto bvpOp = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            nonlinearColumnOperator->getOperator( id ) );
        auto mechOp = std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            bvpOp->getVolumeOperator() );
        mechOp->setVector( AMP::Operator::Mechanics::TEMPERATURE, finalTemperatureVec );
    } // end for id
}

void helperBuildColumnSolverForPelletMechanics(
    std::shared_ptr<AMP::Database> columnSolver_db,
    std::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator,
    std::shared_ptr<AMP::Solver::ColumnSolver> &columnSolver )
{
    auto ikspSolver_db = columnSolver_db->getDatabase( "KrylovSolver" );
    auto mlSolver_db   = ikspSolver_db->getDatabase( "MLSolver" );
    auto columnSolverParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( columnSolver_db );
    columnSolverParams->d_pOperator = linearColumnOperator;
    columnSolver.reset( new AMP::Solver::ColumnSolver( columnSolverParams ) );
    for ( unsigned int id = 0; id < linearColumnOperator->getNumberOfOperators(); id++ ) {
        auto currOp = linearColumnOperator->getOperator( id );
        auto mlSolverParams =
            std::make_shared<AMP::Solver::SolverStrategyParameters>( mlSolver_db );
        mlSolverParams->d_pOperator = currOp;
        auto mlSolver = std::make_shared<AMP::Solver::TrilinosMLSolver>( mlSolverParams );

#ifdef USE_EXT_PETSC
        auto ikspSolverParams =
            std::make_shared<AMP::Solver::PetscKrylovSolverParameters>( ikspSolver_db );
        ikspSolverParams->d_pOperator       = currOp;
        ikspSolverParams->d_comm            = ( currOp->getMesh() )->getComm();
        ikspSolverParams->d_pPreconditioner = mlSolver;
        auto ikspSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( ikspSolverParams );
        columnSolver->append( ikspSolver );
#else
        AMP_ERROR( "petsc required" );
#endif

    } // end for id
}

void helperBuildStackSolverForPelletMechanics(
    std::shared_ptr<AMP::Database> pelletStackSolver_db,
    std::shared_ptr<AMP::Operator::PelletStackOperator> pelletStackOp,
    std::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator,
    std::shared_ptr<AMP::Solver::SolverStrategy> &pelletStackSolver )
{
    auto columnSolver_db = pelletStackSolver_db->getDatabase( "ColumnSolver" );
    std::shared_ptr<AMP::Solver::ColumnSolver> columnSolver;
    helperBuildColumnSolverForPelletMechanics(
        columnSolver_db, linearColumnOperator, columnSolver );
    auto pelletStackSolverParams =
        std::make_shared<AMP::Solver::PelletStackMechanicsSolverParameters>( pelletStackSolver_db );
    pelletStackSolverParams->d_columnSolver = columnSolver;
    pelletStackSolverParams->d_pOperator    = pelletStackOp;
    pelletStackSolver.reset(
        new AMP::Solver::PelletStackMechanicsSolver( pelletStackSolverParams ) );
}

void helperResetNonlinearOperatorForPelletMechanics(
    std::shared_ptr<AMP::Operator::CoupledOperator> coupledOp )
{
    auto nonlinearColumnOperator =
        std::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    auto tmp_db = std::make_shared<AMP::Database>( "Dummy" );
    auto tmpParams =
        std::make_shared<AMP::Operator::MechanicsNonlinearFEOperatorParameters>( tmp_db );
    for ( unsigned int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++ ) {
        auto nonlinOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            nonlinearColumnOperator->getOperator( id ) );
        ( nonlinOperator->getVolumeOperator() )->reset( tmpParams );
    } // end for id
}
