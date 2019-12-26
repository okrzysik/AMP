#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/libmesh/libmeshMesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"
#include "AMP/solvers/petsc/PetscSNESSolver.h"
#include "AMP/solvers/petsc/PetscSNESSolverParameters.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/ReadTestMesh.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/VectorBuilder.h"
#include "libmesh/mesh_communication.h"

#include <iostream>
#include <string>

static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "log_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    auto mesh_file = input_db->getString( "mesh_file" );
        libMesh::Parallel::Communicator comm( AMP_COMM_WORLD );
    auto mesh      = std::make_shared<libMesh::Mesh>( comm, 3 );
    AMP::readTestMesh( mesh_file, mesh );
    libMesh::MeshCommunication().broadcast( *( mesh.get() ) );
    mesh->prepare_for_use( false );

    auto meshAdapter = std::make_shared<AMP::Mesh::libmeshMesh>( mesh, "TestMesh" );

    AMP_INSIST( input_db->keyExists( "NumberOfLoadingSteps" ),
                "Key ''NumberOfLoadingSteps'' is missing!" );
    int NumberOfLoadingSteps = input_db->getScalar<int>( "NumberOfLoadingSteps" );

    AMP_INSIST( input_db->keyExists( "OutputFileName" ), "Key ''OutputFileName'' is missing!" );
    auto outFileName = input_db->getString( "OutputFileName" );

    auto fp = fopen( outFileName.c_str(), "w" );
    fprintf( fp, "clc; \n clear; \n A = zeros(24, 24); \n \n" );

    // Create a nonlinear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "NonlinearMechanicsOperator" ), "key missing!" );
    auto nonlinearMechanicsBVPoperator =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "NonlinearMechanicsOperator", input_db ) );
    auto nonlinearMechanicsVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsBVPoperator->getVolumeOperator() );
    auto mechanicsMaterialModel = nonlinearMechanicsVolumeOperator->getMaterialModel();

    // Create a Linear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "LinearMechanicsOperator" ), "key missing!" );
    // auto linearMechanicsDatabase = input_db->getDatabase("LinearMechanicsOperator");
    auto linearMechanicsBVPoperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "LinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    // Create the variables
    auto mechanicsNonlinearVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsBVPoperator->getVolumeOperator() );
    AMP::LinearAlgebra::Variable::shared_ptr dispVar =
        mechanicsNonlinearVolumeOperator->getOutputVariable();

    // For RHS (Point Forces)
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    auto dirichletLoadVecOp = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "Load_Boundary", input_db, dummyModel ) );
    dirichletLoadVecOp->setVariable( dispVar );

    auto dofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    // Create the vectors
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto solVec       = AMP::LinearAlgebra::createVector( dofMap, dispVar, true );
    auto rhsVec       = solVec->cloneVector();
    auto resVec       = solVec->cloneVector();
    auto scaledRhsVec = solVec->cloneVector();

    // Initial guess
    solVec->zero();

    // RHS
    rhsVec->zero();
    dirichletLoadVecOp->apply( nullVec, rhsVec );
    nonlinearMechanicsBVPoperator->modifyRHSvector( rhsVec );

    // We need to reset the linear operator before the solve since TrilinosML does
    // the factorization of the matrix during construction and so the matrix must
    // be correct before constructing the TrilinosML object.
    nonlinearMechanicsBVPoperator->apply( solVec, resVec );
    linearMechanicsBVPoperator->reset(
        nonlinearMechanicsBVPoperator->getParameters( "Jacobian", solVec ) );

    double epsilon =
        1.0e-13 * ( ( ( linearMechanicsBVPoperator->getMatrix() )->extractDiagonal() )->L1Norm() );

    auto nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    auto linearSolver_db    = nonlinearSolver_db->getDatabase( "LinearSolver" );

    // ---- first initialize the preconditioner
    auto pcSolver_db    = linearSolver_db->getDatabase( "Preconditioner" );
    auto pcSolverParams = std::make_shared<AMP::Solver::TrilinosMLSolverParameters>( pcSolver_db );
    pcSolverParams->d_pOperator = linearMechanicsBVPoperator;
    auto pcSolver               = std::make_shared<AMP::Solver::TrilinosMLSolver>( pcSolverParams );

    // HACK to prevent a double delete on Petsc Vec
    std::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver;

    // initialize the linear solver
    auto linearSolverParams =
        std::make_shared<AMP::Solver::PetscKrylovSolverParameters>( linearSolver_db );
    linearSolverParams->d_pOperator       = linearMechanicsBVPoperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = pcSolver;
    auto linearSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( linearSolverParams );

    // initialize the nonlinear solver
    auto nonlinearSolverParams =
        std::make_shared<AMP::Solver::PetscSNESSolverParameters>( nonlinearSolver_db );
    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nonlinearMechanicsBVPoperator;
    nonlinearSolverParams->d_pKrylovSolver = linearSolver;
    nonlinearSolverParams->d_pInitialGuess = solVec;
    nonlinearSolver.reset( new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );

    nonlinearSolver->setZeroInitialGuess( false );
    /*
       FILE *fout1;
       fout1 = fopen("Loading_Loop.txt","w");
       */
    // double delta_displacement_shear = 0.00006;
    // double delta_displacement_axial = 0.00006;

    for ( int step = 0; step < NumberOfLoadingSteps; step++ ) {
        AMP::pout << "########################################" << std::endl;
        AMP::pout << "The current loading step is " << ( step + 1 ) << std::endl;

        nonlinearMechanicsBVPoperator->modifyInitialSolutionVector( solVec );

        double scaleValue;
        /*    double No5 = NumberOfLoadingSteps / 5;
              double N2o5 = (2 * NumberOfLoadingSteps) / 5;
              double N3o5 = (3 * NumberOfLoadingSteps) / 5;
              double N4o5 = (4 * NumberOfLoadingSteps) / 5;
              double N = NumberOfLoadingSteps;
              if(step < No5) {
              scaleValue = ((double)step+1.0) / ((double)No5);
              }
              if((step >= No5) && (step < N2o5)) {
              scaleValue = 1.0 - (((double)step + 1.0 - (double)No5) / ((double)No5));
              }
              if((step >= N2o5) && (step < N3o5)) {
              scaleValue = - (((double)step + 1.0 - (double)N2o5) / ((double)No5));
              }
              if((step >= N3o5) && (step < N4o5)) {
              scaleValue = -1.0 + (((double)step + 1.0 - (double)N3o5) / ((double)No5));
              }
              if((step >= N4o5) && (step < N)) {
              scaleValue = (((double)step + 1.0 - (double)N4o5) / ((double)No5));
              }
              fprintf(fout1,"%lf %lf\n",((double)step + 1.0),scaleValue);
              */
        /*    double No2 = NumberOfLoadingSteps / 2;
              double N = NumberOfLoadingSteps;
              if(step < No2) {
              scaleValue = ((double)step+1.0) / ((double)No2);
              }
              if((step >= No2) && (step < N)) {
              scaleValue = 1.0 - (((double)step + 1.0 - (double)No2) / ((double)No2));
              }
              fprintf(fout1,"%lf %lf\n",((double)step + 1.0),scaleValue);
              */
        // double No4 = NumberOfLoadingSteps / 4;
        // double No2 = NumberOfLoadingSteps / 2;
        // double N3o4 = (3 * NumberOfLoadingSteps) / 4;
        // double N = NumberOfLoadingSteps;

        scaleValue = ( (double) step + 1.0 ) / NumberOfLoadingSteps;
        scaledRhsVec->scale( scaleValue, rhsVec );
        AMP::pout << "L2 Norm of RHS at loading step " << ( step + 1 ) << " is "
                  << scaledRhsVec->L2Norm() << std::endl;

        nonlinearMechanicsBVPoperator->residual( scaledRhsVec, solVec, resVec );
        double initialResidualNorm = resVec->L2Norm();
        AMP::pout << "Initial Residual Norm for loading step " << ( step + 1 ) << " is "
                  << initialResidualNorm << std::endl;

        nonlinearSolver->solve( scaledRhsVec, solVec );

        nonlinearMechanicsBVPoperator->residual( scaledRhsVec, solVec, resVec );
        double finalResidualNorm = resVec->L2Norm();
        AMP::pout << "Final Residual Norm for loading step " << ( step + 1 ) << " is "
                  << finalResidualNorm << std::endl;

        if ( finalResidualNorm > ( 1.0e-10 * initialResidualNorm ) ) {
            ut->failure( "Nonlinear solve for current loading step" );
        } else {
            ut->passes( "Nonlinear solve for current loading step" );
        }

        double finalSolNorm = solVec->L2Norm();

        AMP::pout << "Final Solution Norm: " << finalSolNorm << std::endl;

        auto mechUvec = solVec->select( AMP::LinearAlgebra::VS_Stride( 0, 3 ), "U" );
        auto mechVvec = solVec->select( AMP::LinearAlgebra::VS_Stride( 1, 3 ), "V" );
        auto mechWvec = solVec->select( AMP::LinearAlgebra::VS_Stride( 2, 3 ), "W" );

        double finalMaxU = mechUvec->maxNorm();
        double finalMaxV = mechVvec->maxNorm();
        double finalMaxW = mechWvec->maxNorm();

        AMP::pout << "Maximum U displacement: " << finalMaxU << std::endl;
        AMP::pout << "Maximum V displacement: " << finalMaxV << std::endl;
        AMP::pout << "Maximum W displacement: " << finalMaxW << std::endl;

        auto tmp_db = std::make_shared<AMP::Database>( "Dummy" );
        auto tmpParams =
            std::make_shared<AMP::Operator::MechanicsNonlinearFEOperatorParameters>( tmp_db );
        ( nonlinearMechanicsBVPoperator->getVolumeOperator() )->reset( tmpParams );

        /*    if(step < (No4 - 1)) {
              dirichletVectorCorrectionDatabase->putScalar("value_1_0", (((double)(step +
           2))*delta_displacement_axial));
              dirichletVectorCorrectionDatabase->putScalar("value_2_0", 0.0);
              std::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new
              AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
              (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
              fprintf(fout1,"%d %le %le\n",step,(((double)(step +
           2))*delta_displacement_axial),0.0);
              }
              if((step >= (No4 - 1)) && (step < (No2 - 1))) {
              dirichletVectorCorrectionDatabase->putScalar("value_1_0", 0.6);
              dirichletVectorCorrectionDatabase->putScalar("value_2_0", (((double)(step + 2 -
           No4))*delta_displacement_shear));
              std::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new
              AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
              (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
              fprintf(fout1,"%d %le %le\n",step,0.3,(((double)(step + 2 -
           No4))*delta_displacement_shear));
              }
              if((step >= (No2 - 1)) && (step < (N3o4 - 1))) {
              dirichletVectorCorrectionDatabase->putScalar("value_1_0", (0.6 - (((double)(step + 2 -
           No2))*delta_displacement_axial)));
              dirichletVectorCorrectionDatabase->putScalar("value_2_0", 0.6);
              std::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new
              AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
              (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
              fprintf(fout1,"%d %le %le\n",step,(0.3 - (((double)(step + 2 -
           No2))*delta_displacement_axial)),0.3);
              }
              if((step >= (N3o4 - 1)) && (step < N)) {
              dirichletVectorCorrectionDatabase->putScalar("value_1_0", 0.0);
              dirichletVectorCorrectionDatabase->putScalar("value_2_0", (0.6 - (((double)(step + 2 -
           N3o4))*delta_displacement_shear)));
              std::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new
              AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
              (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
              fprintf(fout1,"%d %le %le\n",step,0.0,(0.3 - (((double)(step + 2 -
           N3o4))*delta_displacement_shear)));
              }
              */

        /*    if(step == 0) {
              dirichletVectorCorrectionDatabase->putScalar("value_1_0", 0.3);
              dirichletVectorCorrectionDatabase->putScalar("value_2_0", 0.3);
              std::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new
              AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
              (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
              fprintf(fout1,"%d %le %le\n",step,(((double)(step + 2))*0.003),0.0);
              }
              if(step == 1) {
              dirichletVectorCorrectionDatabase->putScalar("value_1_0", 0.0);
              dirichletVectorCorrectionDatabase->putScalar("value_2_0", 0.3);
              std::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new
              AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
              (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
              fprintf(fout1,"%d %le %le\n",step,0.3,(((double)(step + 2 - No4))*0.003));
              }
              if((step == 2) || (step == 3)) {
              dirichletVectorCorrectionDatabase->putScalar("value_1_0", 0.0);
              dirichletVectorCorrectionDatabase->putScalar("value_2_0", 0.0);
              std::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> bndParams(new
              AMP::Operator::DirichletVectorCorrectionParameters(dirichletVectorCorrectionDatabase));
              (nonlinearMechanicsBVPoperator->getBoundaryOperator())->reset(bndParams);
              fprintf(fout1,"%d %le %le\n",step,(0.3 - (((double)(step + 2 - No2))*0.003)),0.3);
              }
              */
        nonlinearSolver->setZeroInitialGuess( false );

        // std::cout<<solVec<<std::endl;

        auto mechMat = linearMechanicsBVPoperator->getMatrix();

        for ( int i = 0; i < 24; i++ ) {
            std::vector<size_t> matCols;
            std::vector<double> matVals;
            mechMat->getRowByGlobalID( i, matCols, matVals );
            for ( unsigned int j = 0; j < matCols.size(); j++ ) {
                fprintf(
                    fp, "A(%d, %d) = %.15f ; \n", ( i + 1 ), (int) ( matCols[j] + 1 ), matVals[j] );
            } // end for j
            fprintf( fp, "\n" );
        } // end for i

        /*char num1[256];
          sprintf(num1,"%d",step);
          auto number1 = num1;
          auto fname = exeName + "_Stress_Strain_" + number1 + ".txt";
          std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(nonlinearMechanicsBVPoperator->getVolumeOperator())->printStressAndStrain(solVec,
          fname);
          */
    }

    // fclose(fout1);

    // AMP::pout<<solVec<<std::endl;

    AMP::pout << "epsilon = " << epsilon << std::endl;

    // mechanicsNonlinearVolumeOperator->printStressAndStrain(solVec, output_file);

    ut->passes( exeName );
}

int testUpdatedLagrangianMechanics_LinearElasticity_1( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    auto libmeshInit = std::make_shared<AMP::Mesh::initializeLibMesh>( AMP_COMM_WORLD );

    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    // exeNames.push_back("testUpdatedLagrangianMechanics-LinearElasticity-1");
    // exeNames.push_back("testUpdatedLagrangianMechanics-LinearElasticity-1a");
    exeNames.emplace_back( "testUpdatedLagrangianMechanics-LinearElasticity-1b" );
    // exeNames.push_back("testUpdatedLagrangianMechanics-LinearElasticity-1c");

    for ( auto &exeName : exeNames )
        myTest( &ut, exeName );

    ut.report();
    int num_failed = ut.NumFailGlobal();

    libmeshInit.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
