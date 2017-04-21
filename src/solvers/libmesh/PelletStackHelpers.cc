#include "solvers/libmesh/PelletStackHelpers.h"

#include "discretization/simpleDOF_Manager.h"
#include "matrices/MatrixBuilder.h"
#include "operators/CoupledOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/libmesh/PelletStackOperator.h"
#include "operators/map/AsyncMapColumnOperator.h"
#include "operators/map/NodeToNodeMap.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "vectors/VectorBuilder.h"

#ifdef USE_EXT_PETSC
#include "solvers/petsc/PetscKrylovSolver.h"
#endif

#include "solvers/trilinos/ml/TrilinosMLSolver.h"


void helperCreateStackOperatorForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator>
        n2nmaps,
    AMP::shared_ptr<AMP::Database>
        global_input_db,
    AMP::shared_ptr<AMP::Operator::PelletStackOperator> &pelletStackOp )
{
    AMP::shared_ptr<AMP::Database> pelletStackOp_db =
        global_input_db->getDatabase( "PelletStackOperator" );
    AMP::shared_ptr<AMP::Mesh::Mesh> pelletMeshes = manager->Subset( "PelletMeshes" );
    std::vector<AMP::Mesh::MeshID> pelletMeshIDs  = pelletMeshes->getBaseMeshIDs();
    unsigned int totalNumPellets                  = pelletMeshIDs.size();
    pelletStackOp_db->putInteger( "TOTAL_NUMBER_OF_PELLETS", totalNumPellets );
    AMP::shared_ptr<AMP::Operator::PelletStackOperatorParameters> pelletStackOpParams(
        new AMP::Operator::PelletStackOperatorParameters( pelletStackOp_db ) );
    pelletStackOpParams->d_pelletStackComm = pelletMeshes->getComm();
    pelletStackOpParams->d_n2nMaps         = n2nmaps;
    pelletStackOpParams->d_Mesh            = pelletMeshes;
    pelletStackOp.reset( new AMP::Operator::PelletStackOperator( pelletStackOpParams ) );
}

void helperCreateColumnOperatorsForPelletMechanics(
    std::vector<unsigned int> localPelletIds,
    std::vector<AMP::Mesh::Mesh::shared_ptr>
        localMeshes,
    AMP::shared_ptr<AMP::Database>
        global_input_db,
    AMP::shared_ptr<AMP::Operator::ColumnOperator> &nonlinearColumnOperator,
    AMP::shared_ptr<AMP::Operator::ColumnOperator> &linearColumnOperator )
{
    AMP::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
    nonlinearColumnOperator.reset( new AMP::Operator::ColumnOperator( emptyParams ) );
    linearColumnOperator.reset( new AMP::Operator::ColumnOperator( emptyParams ) );
    for ( unsigned int id = 0; id < localPelletIds.size(); id++ ) {
        std::string prefix = "";
        if ( localPelletIds[id] == 0 ) {
            prefix = "Bottom";
        }

        AMP::Mesh::Mesh::shared_ptr meshAdapter = localMeshes[id];

        AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> mechModel;
        AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinOperator =
            AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter,
                    prefix + "PelletMechanicsNonlinearBVPOperator",
                    global_input_db,
                    mechModel ) );
        nonlinearColumnOperator->append( nonlinOperator );

        AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linOperator =
            AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter,
                    prefix + "PelletMechanicsLinearBVPOperator",
                    global_input_db,
                    mechModel ) );
        linearColumnOperator->append( linOperator );
    } // end for id
}

void helperCreateCoupledOperatorForPelletMechanics(
    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps,
    AMP::shared_ptr<AMP::Operator::ColumnOperator>
        nonlinearColumnOperator,
    AMP::shared_ptr<AMP::Operator::CoupledOperator> &coupledOp )
{
    AMP::shared_ptr<AMP::Database> emptyDb;
    AMP::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledOpParams(
        new AMP::Operator::CoupledOperatorParameters( emptyDb ) );
    coupledOpParams->d_MapOperator = n2nmaps;
    coupledOpParams->d_BVPOperator = nonlinearColumnOperator;
    coupledOp.reset( new AMP::Operator::CoupledOperator( coupledOpParams ) );
}

void helperSetFrozenVectorForMapsForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager, AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledOp )
{
    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps =
        AMP::dynamic_pointer_cast<AMP::Operator::AsyncMapColumnOperator>(
            coupledOp->getOperator( 2 ) );
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    AMP::LinearAlgebra::Variable::shared_ptr dispVar = nonlinearColumnOperator->getOutputVariable();
    AMP::Discretization::DOFManager::shared_ptr nodal3VectorDOF =
        AMP::Discretization::simpleDOFManager::create( manager, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    AMP::LinearAlgebra::Vector::shared_ptr dirichletValues =
        AMP::LinearAlgebra::createVector( nodal3VectorDOF, dispVar, true );
    if ( n2nmaps ) {
        n2nmaps->setVector( dirichletValues );
    }
    for ( unsigned int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++ ) {
        AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletOp =
            AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
                AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
                    nonlinearColumnOperator->getOperator( id ) )
                    ->getBoundaryOperator() );
        dirichletOp->setDirichletValues( dirichletValues );
    } // end for id
}

void helperCreateAllOperatorsForPelletMechanics(
    AMP::Mesh::Mesh::shared_ptr manager,
    AMP::AMP_MPI,
    AMP::shared_ptr<AMP::Database>
        global_input_db,
    AMP::shared_ptr<AMP::Operator::CoupledOperator> &coupledOp,
    AMP::shared_ptr<AMP::Operator::ColumnOperator> &linearColumnOperator,
    AMP::shared_ptr<AMP::Operator::PelletStackOperator> &pelletStackOp )
{
    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps;
    if ( global_input_db->keyExists( "MechanicsNodeToNodeMaps" ) ) {
        AMP::shared_ptr<AMP::Database> map_db =
            global_input_db->getDatabase( "MechanicsNodeToNodeMaps" );
        n2nmaps = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap>(
            manager, map_db );
    }
    helperCreateStackOperatorForPelletMechanics( manager, n2nmaps, global_input_db, pelletStackOp );
    std::vector<unsigned int> localPelletIds             = pelletStackOp->getLocalPelletIds();
    std::vector<AMP::Mesh::Mesh::shared_ptr> localMeshes = pelletStackOp->getLocalMeshes();
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator;
    helperCreateColumnOperatorsForPelletMechanics( localPelletIds,
                                                   localMeshes,
                                                   global_input_db,
                                                   nonlinearColumnOperator,
                                                   linearColumnOperator );
    helperCreateCoupledOperatorForPelletMechanics( n2nmaps, nonlinearColumnOperator, coupledOp );
    helperSetFrozenVectorForMapsForPelletMechanics( manager, coupledOp );
}

void helperCreateVectorsForPelletMechanics( AMP::Mesh::Mesh::shared_ptr manager,
                                            AMP::shared_ptr<AMP::Operator::CoupledOperator>
                                                coupledOp,
                                            AMP::LinearAlgebra::Vector::shared_ptr &solVec,
                                            AMP::LinearAlgebra::Vector::shared_ptr &rhsVec,
                                            AMP::LinearAlgebra::Vector::shared_ptr &scaledRhsVec )
{
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    AMP::LinearAlgebra::Variable::shared_ptr dispVar = nonlinearColumnOperator->getOutputVariable();
    AMP::Discretization::DOFManager::shared_ptr nodal3VectorDOF =
        AMP::Discretization::simpleDOFManager::create( manager, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    solVec       = AMP::LinearAlgebra::createVector( nodal3VectorDOF, dispVar, true );
    rhsVec       = AMP::LinearAlgebra::createVector( nodal3VectorDOF, dispVar, true );
    scaledRhsVec = rhsVec->cloneVector();
}

void helperBuildPointLoadRHSForPelletMechanics( AMP::shared_ptr<AMP::Database> global_input_db,
                                                AMP::shared_ptr<AMP::Operator::CoupledOperator>
                                                    coupledOp,
                                                AMP::LinearAlgebra::Vector::shared_ptr rhsVec )
{
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    rhsVec->zero();
    for ( unsigned int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++ ) {
        AMP::Operator::Operator::shared_ptr currOp = nonlinearColumnOperator->getOperator( id );
        AMP::Mesh::Mesh::shared_ptr meshAdapter    = currOp->getMesh();
        AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
        AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOp =
            AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter, "PointLoad", global_input_db, dummyModel ) );
        loadOp->setVariable( currOp->getOutputVariable() );
        loadOp->apply( nullVec, rhsVec );
    } // end for id
}

void helperApplyBoundaryCorrectionsForPelletMechanics(
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr solVec,
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec )
{
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    for ( unsigned int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++ ) {
        AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinOperator =
            AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
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
    AMP::LinearAlgebra::Variable::shared_ptr tempVar(
        new AMP::LinearAlgebra::Variable( "temperature" ) );
    AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF =
        AMP::Discretization::simpleDOFManager::create( manager, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    initialTemperatureVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, tempVar, true );
    finalTemperatureVec   = initialTemperatureVec->cloneVector();
}

void helperSetReferenceTemperatureForPelletMechanics(
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr initialTemperatureVec )
{
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    for ( unsigned int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++ ) {
        AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> bvpOp =
            AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
                nonlinearColumnOperator->getOperator( id ) );
        AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechOp =
            AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
                bvpOp->getVolumeOperator() );
        mechOp->setReferenceTemperature( initialTemperatureVec );
    } // end for id
}

void helperSetFinalTemperatureForPelletMechanics(
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr finalTemperatureVec )
{
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    for ( unsigned int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++ ) {
        AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> bvpOp =
            AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
                nonlinearColumnOperator->getOperator( id ) );
        AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechOp =
            AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
                bvpOp->getVolumeOperator() );
        mechOp->setVector( AMP::Operator::Mechanics::TEMPERATURE, finalTemperatureVec );
    } // end for id
}

void helperBuildColumnSolverForPelletMechanics(
    AMP::shared_ptr<AMP::Database> columnSolver_db,
    AMP::shared_ptr<AMP::Operator::ColumnOperator>
        linearColumnOperator,
    AMP::shared_ptr<AMP::Solver::ColumnSolver> &columnSolver )
{
    AMP::shared_ptr<AMP::Database> ikspSolver_db = columnSolver_db->getDatabase( "KrylovSolver" );
    AMP::shared_ptr<AMP::Database> mlSolver_db   = ikspSolver_db->getDatabase( "MLSolver" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> columnSolverParams(
        new AMP::Solver::SolverStrategyParameters( columnSolver_db ) );
    columnSolverParams->d_pOperator = linearColumnOperator;
    columnSolver.reset( new AMP::Solver::ColumnSolver( columnSolverParams ) );
    for ( unsigned int id = 0; id < linearColumnOperator->getNumberOfOperators(); id++ ) {
        AMP::Operator::Operator::shared_ptr currOp = linearColumnOperator->getOperator( id );
        AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams(
            new AMP::Solver::SolverStrategyParameters( mlSolver_db ) );
        mlSolverParams->d_pOperator = currOp;
        AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(
            new AMP::Solver::TrilinosMLSolver( mlSolverParams ) );

#ifdef USE_EXT_PETSC
        AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> ikspSolverParams(
            new AMP::Solver::PetscKrylovSolverParameters( ikspSolver_db ) );
        ikspSolverParams->d_pOperator       = currOp;
        ikspSolverParams->d_comm            = ( currOp->getMesh() )->getComm();
        ikspSolverParams->d_pPreconditioner = mlSolver;
        AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> ikspSolver(
            new AMP::Solver::PetscKrylovSolver( ikspSolverParams ) );
        columnSolver->append( ikspSolver );
#else
        AMP_ERROR( "petsc required" );
#endif

    } // end for id
}

void helperBuildStackSolverForPelletMechanics(
    AMP::shared_ptr<AMP::Database> pelletStackSolver_db,
    AMP::shared_ptr<AMP::Operator::PelletStackOperator>
        pelletStackOp,
    AMP::shared_ptr<AMP::Operator::ColumnOperator>
        linearColumnOperator,
    AMP::shared_ptr<AMP::Solver::SolverStrategy> &pelletStackSolver )
{
    AMP::shared_ptr<AMP::Database> columnSolver_db =
        pelletStackSolver_db->getDatabase( "ColumnSolver" );
    AMP::shared_ptr<AMP::Solver::ColumnSolver> columnSolver;
    helperBuildColumnSolverForPelletMechanics(
        columnSolver_db, linearColumnOperator, columnSolver );
    AMP::shared_ptr<AMP::Solver::PelletStackMechanicsSolverParameters> pelletStackSolverParams(
        new AMP::Solver::PelletStackMechanicsSolverParameters( pelletStackSolver_db ) );
    pelletStackSolverParams->d_columnSolver = columnSolver;
    pelletStackSolverParams->d_pOperator    = pelletStackOp;
    pelletStackSolver.reset(
        new AMP::Solver::PelletStackMechanicsSolver( pelletStackSolverParams ) );
}

void helperResetNonlinearOperatorForPelletMechanics(
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledOp )
{
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( coupledOp->getOperator( 3 ) );
    AMP::shared_ptr<AMP::Database> tmp_db( new AMP::InputDatabase( "Dummy" ) );
    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> tmpParams(
        new AMP::Operator::MechanicsNonlinearFEOperatorParameters( tmp_db ) );

    for ( unsigned int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++ ) {
        AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinOperator =
            AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
                nonlinearColumnOperator->getOperator( id ) );
        ( nonlinOperator->getVolumeOperator() )->reset( tmpParams );
    } // end for id
}
