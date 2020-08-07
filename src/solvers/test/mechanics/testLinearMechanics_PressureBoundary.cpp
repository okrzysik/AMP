#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/PressureBoundaryOperator.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/ReadTestMesh.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"

#include <fstream>
#include <iostream>
#include <string>


static void linearElasticTest( AMP::UnitTest *ut, std::string exeName, int exampleNum )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName + ".txt";

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

#ifdef USE_EXT_SILO
    // Create the silo writer and register the data
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
#endif


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    std::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "MechanicsBVPOperator", input_db, elementPhysicsModel ) );

    AMP::LinearAlgebra::Variable::shared_ptr dispVar = bvpOperator->getOutputVariable();

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    std::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletVecOp =
        std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "Load_Boundary", input_db, dummyModel ) );
    // This has an in-place apply. So, it has an empty input variable and
    // the output variable is the same as what it is operating on.
    dirichletVecOp->setVariable( dispVar );

    // Pressure RHS
    std::shared_ptr<AMP::Operator::PressureBoundaryOperator> pressureLoadVecOp =
        std::dynamic_pointer_cast<AMP::Operator::PressureBoundaryOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "Pressure_Boundary", input_db, dummyModel ) );
    // This has an in-place apply. So, it has an empty input variable and
    // the output variable is the same as what it is operating on.

    AMP::Discretization::DOFManager::shared_ptr dofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Vector::shared_ptr mechSolVec =
        AMP::LinearAlgebra::createVector( dofMap, dispVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr mechRhsVec      = mechSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr mechResVec      = mechSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr mechPressureVec = mechSolVec->cloneVector();

    mechSolVec->setToScalar( 0.0, mechSolVec );
    mechRhsVec->setToScalar( 0.0, mechRhsVec );
    mechResVec->setToScalar( 0.0, mechResVec );
    mechPressureVec->setToScalar( 0.0, mechPressureVec );

    dirichletVecOp->apply( nullVec, mechRhsVec );

    double rhsNorm = mechRhsVec->L2Norm( mechRhsVec );
    AMP::pout << "RHS Norm after Dirichlet Apply: " << rhsNorm << std::endl;

    double pressNorm = mechPressureVec->L2Norm( mechPressureVec );
    AMP::pout << "Pressure Norm before Apply: " << pressNorm << std::endl;

    // Applying the pressure load
    pressureLoadVecOp->addRHScorrection( mechPressureVec );

    pressNorm = mechPressureVec->L2Norm( mechPressureVec );
    AMP::pout << "Pressure Norm after Apply: " << pressNorm << std::endl;

    mechRhsVec->add( mechRhsVec, mechPressureVec, mechRhsVec );

    rhsNorm = mechRhsVec->L2Norm( mechRhsVec );
    AMP::pout << "Total RHS Norm: " << rhsNorm << std::endl;

    double initSolNorm = mechSolVec->L2Norm( mechSolVec );

    AMP::pout << "Initial Solution Norm: " << initSolNorm << std::endl;

    bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

    double initResidualNorm = mechResVec->L2Norm( mechResVec );

    AMP::pout << "Initial Residual Norm: " << initResidualNorm << std::endl;

    std::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase( "LinearSolver" );

    // ---- first initialize the preconditioner
    std::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase( "Preconditioner" );
    std::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(
        new AMP::Solver::TrilinosMLSolverParameters( pcSolver_db ) );
    pcSolverParams->d_pOperator = bvpOperator;
    std::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(
        new AMP::Solver::TrilinosMLSolver( pcSolverParams ) );

    // initialize the linear solver
    std::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
        new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
    linearSolverParams->d_pOperator       = bvpOperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = pcSolver;
    std::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
        new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );

    linearSolver->setZeroInitialGuess( false );

    linearSolver->solve( mechRhsVec, mechSolVec );

    double finalSolNorm = mechSolVec->L2Norm( mechSolVec );

    AMP::pout << "Final Solution Norm: " << finalSolNorm << std::endl;

    std::string fname = exeName + "_StressAndStrain.txt";

    ( std::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(
          bvpOperator->getVolumeOperator() ) )
        ->printStressAndStrain( mechSolVec, fname );

    AMP::LinearAlgebra::Vector::shared_ptr mechUvec =
        mechSolVec->select( AMP::LinearAlgebra::VS_Stride( 0, 3 ), "U" );
    AMP::LinearAlgebra::Vector::shared_ptr mechVvec =
        mechSolVec->select( AMP::LinearAlgebra::VS_Stride( 1, 3 ), "V" );
    AMP::LinearAlgebra::Vector::shared_ptr mechWvec =
        mechSolVec->select( AMP::LinearAlgebra::VS_Stride( 2, 3 ), "W" );

    double finalMaxU = mechUvec->maxNorm( mechUvec );
    double finalMaxV = mechVvec->maxNorm( mechVvec );
    double finalMaxW = mechWvec->maxNorm( mechWvec );

    AMP::pout << "Maximum U displacement: " << finalMaxU << std::endl;
    AMP::pout << "Maximum V displacement: " << finalMaxV << std::endl;
    AMP::pout << "Maximum W displacement: " << finalMaxW << std::endl;

    bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

    double finalResidualNorm = mechResVec->L2Norm( mechResVec );

    AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm > ( 1e-10 * initResidualNorm ) ) {
        ut->failure( exeName );
    } else {
        ut->passes( exeName );
    }

    auto diag      = ( bvpOperator->getMatrix() )->extractDiagonal();
    double epsilon = 1.0e-13 * diag->L1Norm( diag );
    AMP::pout << "epsilon = " << epsilon << std::endl;

#ifdef USE_EXT_SILO
    siloWriter->registerVector( mechSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    char outFileName1[256];
    sprintf( outFileName1, "undeformedBeam_%d", exampleNum );
    siloWriter->writeFile( outFileName1, 1 );
    meshAdapter->displaceMesh( mechSolVec );
    char outFileName2[256];
    sprintf( outFileName2, "deformedBeam_%d", exampleNum );
    siloWriter->writeFile( outFileName2, 1 );
#endif
}

int testLinearMechanics_PressureBoundary( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;

    if ( argc == 1 ) {
        exeNames.emplace_back( "testLinearMechanics-PressureBoundary-1" );
        exeNames.emplace_back( "testLinearMechanics-PressureBoundary-HaldenPellet" );
    } else {
        for ( int i = 1; i < argc; ++i ) {
            char inpName[100];
            sprintf( inpName, "testLinearMechanics-PressureBoundary-%s", argv[i] );
            exeNames.emplace_back( inpName );
        } // end for i
    }

    for ( size_t i = 0; i < exeNames.size(); ++i ) {
        linearElasticTest( &ut, exeNames[i], i );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
