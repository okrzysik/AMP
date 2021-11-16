#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscSNESSolver.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"

#include <iostream>
#include <memory>
#include <string>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

#ifdef USE_EXT_SILO
    // Create the silo writer and register the data
    auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
#endif

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db    = input_db->getDatabase( "Mesh" );
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    // create a nonlinear BVP operator for nonlinear mechanics
    AMP_INSIST( input_db->keyExists( "testNonlinearMechanicsOperator" ), "key missing!" );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel;
    auto nonlinearMechanicsOperator =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    // create a nonlinear BVP operator for nonlinear thermal diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearThermalOperator" ), "key missing!" );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
    auto nonlinearThermalOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testNonlinearThermalOperator", input_db, thermalTransportModel ) );

    // create a column operator object for nonlinear thermomechanics
    auto nonlinearThermoMechanicsOperator = std::make_shared<AMP::Operator::ColumnOperator>();
    nonlinearThermoMechanicsOperator->append( nonlinearMechanicsOperator );
    nonlinearThermoMechanicsOperator->append( nonlinearThermalOperator );

    // initialize the input multi-variable
    auto mechanicsVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsOperator->getVolumeOperator() );

    // initialize the output multi-variable
    auto displacementVar = nonlinearMechanicsOperator->getOutputVariable();
    auto temperatureVar  = nonlinearThermalOperator->getOutputVariable();

    auto vectorDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    std::shared_ptr<AMP::Discretization::DOFManager> scalarDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    // create solution, rhs, and residual vectors
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto displacementVec = AMP::LinearAlgebra::createVector( vectorDofMap, displacementVar, true );
    auto temperatureVec  = AMP::LinearAlgebra::createVector( scalarDofMap, temperatureVar, true );
    auto solVec          = AMP::LinearAlgebra::MultiVector::create( "multiVector", globalComm );
    auto multiVec        = std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>( solVec );
    multiVec->addVector( displacementVec );
    multiVec->addVector( temperatureVec );

    auto rhsVec = solVec->cloneVector();
    auto resVec = solVec->cloneVector();

#ifdef USE_EXT_SILO
    siloWriter->registerVector(
        displacementVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "MechanicsSolution" );
    siloWriter->registerVector(
        temperatureVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "ThermalSolution" );
#endif

    auto referenceTemperatureVec = temperatureVec->cloneVector();
    referenceTemperatureVec->setToScalar( 300.0 );
    mechanicsVolumeOperator->setReferenceTemperature( referenceTemperatureVec );

    // now construct the linear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "testLinearMechanicsOperator" ), "key missing!" );
    auto linearMechanicsOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testLinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists( "testLinearThermalOperator" ), "key missing!" );
    auto linearThermalOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testLinearThermalOperator", input_db, thermalTransportModel ) );

    // create a column operator object for linear thermomechanics
    auto linearThermoMechanicsOperator = std::make_shared<AMP::Operator::ColumnOperator>();
    linearThermoMechanicsOperator->append( linearMechanicsOperator );
    linearThermoMechanicsOperator->append( linearThermalOperator );

    // Initial-Guess for mechanics
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyMechanicsModel;
    auto dirichletDispInVecOp = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "MechanicsInitialGuess", input_db, dummyMechanicsModel ) );
    dirichletDispInVecOp->setVariable( displacementVar );

    // Initial-Guess for thermal
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyThermalModel;
    auto dirichletThermalInVecOp =
        std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "ThermalInitialGuess", input_db, dummyThermalModel ) );
    dirichletThermalInVecOp->setVariable( temperatureVar );

    // Random initial guess
    solVec->setToScalar( 0.0 );
    const double referenceTemperature = 301.0;
    temperatureVec->addScalar( *temperatureVec, referenceTemperature );

    // Initial guess for mechanics must satisfy the displacement boundary conditions
    dirichletDispInVecOp->apply( nullVec, solVec );
    // Initial guess for thermal must satisfy the thermal Dirichlet boundary conditions
    dirichletThermalInVecOp->apply( nullVec, solVec );

    // We need to reset the linear operator before the solve since TrilinosML does
    // the factorization of the matrix during construction and so the matrix must
    // be correct before constructing the TrilinosML object.
    // The thermal operator does not expect an apply to be called before calling
    // getJacobianParams and so it need not be called. So, any of the following
    // apply calls will work:
    // mechanicsVolumeOperator->apply(nullVec, solVec, resVec, 1.0, 0.0);
    // nonlinearMechanicsOperator->apply(nullVec, solVec, resVec, 1.0, 0.0);
    nonlinearThermoMechanicsOperator->apply( solVec, resVec );
    linearThermoMechanicsOperator->reset(
        nonlinearThermoMechanicsOperator->getParameters( "Jacobian", solVec ) );

    rhsVec->setToScalar( 0.0 );

    auto nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    auto linearSolver_db    = nonlinearSolver_db->getDatabase( "LinearSolver" );

    // initialize the nonlinear solver
    auto nonlinearSolverParams =
        std::make_shared<AMP::Solver::PetscSNESSolverParameters>( nonlinearSolver_db );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nonlinearThermoMechanicsOperator;
    nonlinearSolverParams->d_pInitialGuess = solVec;
    auto nonlinearSolver = std::make_shared<AMP::Solver::PetscSNESSolver>( nonlinearSolverParams );

    // initialize the column preconditioner which is a diagonal block preconditioner
    auto columnPreconditioner_db = linearSolver_db->getDatabase( "Preconditioner" );
    auto columnPreconditionerParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( columnPreconditioner_db );
    columnPreconditionerParams->d_pOperator = linearThermoMechanicsOperator;
    auto columnPreconditioner =
        std::make_shared<AMP::Solver::ColumnSolver>( columnPreconditionerParams );

    auto mechanicsPreconditioner_db =
        columnPreconditioner_db->getDatabase( "mechanicsPreconditioner" );
    auto mechanicsPreconditionerParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( mechanicsPreconditioner_db );
    mechanicsPreconditionerParams->d_pOperator = linearMechanicsOperator;
    auto linearMechanicsPreconditioner =
        std::make_shared<AMP::Solver::TrilinosMLSolver>( mechanicsPreconditionerParams );

    auto thermalPreconditioner_db = columnPreconditioner_db->getDatabase( "thermalPreconditioner" );
    auto thermalPreconditionerParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( thermalPreconditioner_db );
    thermalPreconditionerParams->d_pOperator = linearThermalOperator;
    auto linearThermalPreconditioner =
        std::make_shared<AMP::Solver::TrilinosMLSolver>( thermalPreconditionerParams );

    columnPreconditioner->append( linearMechanicsPreconditioner );
    columnPreconditioner->append( linearThermalPreconditioner );

    // register the preconditioner with the Jacobian free Krylov solver
    auto linearSolver = nonlinearSolver->getKrylovSolver();

    linearSolver->setPreconditioner( columnPreconditioner );

    nonlinearThermoMechanicsOperator->residual( rhsVec, solVec, resVec );
    double initialResidualNorm = static_cast<double>( resVec->L2Norm() );

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;

    nonlinearSolver->setZeroInitialGuess( false );

    nonlinearSolver->apply( rhsVec, solVec );

    nonlinearThermoMechanicsOperator->residual( rhsVec, solVec, resVec );

    double finalResidualNorm = static_cast<double>( resVec->L2Norm() );
    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    auto mechUvec = displacementVec->select( AMP::LinearAlgebra::VS_Stride( 0, 3 ), "U" );
    auto mechVvec = displacementVec->select( AMP::LinearAlgebra::VS_Stride( 1, 3 ), "V" );
    auto mechWvec = displacementVec->select( AMP::LinearAlgebra::VS_Stride( 2, 3 ), "W" );

    double finalMaxU = static_cast<double>( mechUvec->maxNorm() );
    double finalMaxV = static_cast<double>( mechVvec->maxNorm() );
    double finalMaxW = static_cast<double>( mechWvec->maxNorm() );

    AMP::pout << "Maximum U displacement: " << finalMaxU << std::endl;
    AMP::pout << "Maximum V displacement: " << finalMaxV << std::endl;
    AMP::pout << "Maximum W displacement: " << finalMaxW << std::endl;

#ifdef USE_EXT_SILO
    siloWriter->writeFile( exeName, 1 );
#endif

    if ( finalResidualNorm > initialResidualNorm * 1.0e-10 + 1.0e-05 ) {
        ut->failure( "Error" );
    } else {
        ut->passes( "PetscSNES Solver successfully solves a nonlinear thermo-mechanics equation "
                    "with JFNK, FGMRES for "
                    "Krylov, block diagonal preconditioning with ML solvers" );
    }
    ut->passes( exeName );
}


int testPetscSNESSolver_NonlinearThermoMechanics_1( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testPetscSNESSolver-NonlinearThermoMechanics-1" );
    exeNames.emplace_back( "testPetscSNESSolver-NonlinearThermoMechanics-1a" );
    exeNames.emplace_back( "testPetscSNESSolver-NonlinearThermoMechanics-1b" );
    exeNames.emplace_back( "testPetscSNESSolver-NonlinearThermoMechanics-1c" );

    for ( auto &exeName : exeNames )
        myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
