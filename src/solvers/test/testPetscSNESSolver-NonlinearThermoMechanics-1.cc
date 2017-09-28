
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "discretization/simpleDOF_Manager.h"
#include "utils/Writer.h"
#include "vectors/VectorBuilder.h"

#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"

#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "solvers/ColumnSolver.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscSNESSolver.h"
#include "solvers/trilinos/ml/TrilinosMLSolver.h"

void myTest( AMP::UnitTest *ut, const std::string& exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

#ifdef USE_EXT_SILO
    // Create the silo writer and register the data
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
#endif

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear mechanics
    AMP_INSIST( input_db->keyExists( "testNonlinearMechanicsOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear thermal diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearThermalOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearThermalOperator", input_db, thermalTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a column operator object for nonlinear thermomechanics
    AMP::shared_ptr<AMP::Operator::OperatorParameters> params;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearThermoMechanicsOperator(
        new AMP::Operator::ColumnOperator( params ) );
    nonlinearThermoMechanicsOperator->append( nonlinearMechanicsOperator );
    nonlinearThermoMechanicsOperator->append( nonlinearThermalOperator );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the input multi-variable
    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechanicsVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsOperator->getVolumeOperator() );

    // initialize the output multi-variable
    AMP::LinearAlgebra::Variable::shared_ptr displacementVar =
        nonlinearMechanicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr temperatureVar =
        nonlinearThermalOperator->getOutputVariable();

    AMP::Discretization::DOFManager::shared_ptr vectorDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    AMP::Discretization::DOFManager::shared_ptr scalarDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    // create solution, rhs, and residual vectors
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Vector::shared_ptr displacementVec =
        AMP::LinearAlgebra::createVector( vectorDofMap, displacementVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr temperatureVec =
        AMP::LinearAlgebra::createVector( scalarDofMap, temperatureVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::MultiVector::create( "multiVector", globalComm );
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> multiVec =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>( solVec );
    multiVec->addVector( displacementVec );
    multiVec->addVector( temperatureVec );

    AMP::LinearAlgebra::Vector::shared_ptr rhsVec = solVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr resVec = solVec->cloneVector();

    //----------------------------------------------------------------------------------------------------------------------------------------------//

#ifdef USE_EXT_SILO
    siloWriter->registerVector(
        displacementVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "MechanicsSolution" );
    siloWriter->registerVector(
        temperatureVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "ThermalSolution" );
#endif

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    AMP::LinearAlgebra::Vector::shared_ptr referenceTemperatureVec = temperatureVec->cloneVector();
    referenceTemperatureVec->setToScalar( 300.0 );
    mechanicsVolumeOperator->setReferenceTemperature( referenceTemperatureVec );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "testLinearMechanicsOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists( "testLinearThermalOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearThermalOperator", input_db, thermalTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a column operator object for linear thermomechanics
    AMP::shared_ptr<AMP::Operator::ColumnOperator> linearThermoMechanicsOperator(
        new AMP::Operator::ColumnOperator( params ) );
    linearThermoMechanicsOperator->append( linearMechanicsOperator );
    linearThermoMechanicsOperator->append( linearThermalOperator );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Initial-Guess for mechanics
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyMechanicsModel;
    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletDispInVecOp =
        AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "MechanicsInitialGuess", input_db, dummyMechanicsModel ) );
    dirichletDispInVecOp->setVariable( displacementVar );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Initial-Guess for thermal
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyThermalModel;
    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletThermalInVecOp =
        AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "ThermalInitialGuess", input_db, dummyThermalModel ) );
    dirichletThermalInVecOp->setVariable( temperatureVar );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Random initial guess
    solVec->setToScalar( 0.0 );
    const double referenceTemperature = 301.0;
    temperatureVec->addScalar( temperatureVec, referenceTemperature );

    // Initial guess for mechanics must satisfy the displacement boundary conditions
    dirichletDispInVecOp->apply( nullVec, solVec );
    // Initial guess for thermal must satisfy the thermal Dirichlet boundary conditions
    dirichletThermalInVecOp->apply( nullVec, solVec );

    //----------------------------------------------------------------------------------------------------------------------------------------------//

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
    //----------------------------------------------------------------------------------------------------------------------------------------------/

    rhsVec->setToScalar( 0.0 );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );

    AMP::shared_ptr<AMP::Database> linearSolver_db =
        nonlinearSolver_db->getDatabase( "LinearSolver" );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the nonlinear solver
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db ) );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm      = globalComm;
    nonlinearSolverParams->d_pOperator = nonlinearThermoMechanicsOperator;

    nonlinearSolverParams->d_pInitialGuess = solVec;

    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(
        new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the column preconditioner which is a diagonal block preconditioner
    AMP::shared_ptr<AMP::Database> columnPreconditioner_db =
        linearSolver_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( columnPreconditioner_db ) );
    columnPreconditionerParams->d_pOperator = linearThermoMechanicsOperator;
    AMP::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(
        new AMP::Solver::ColumnSolver( columnPreconditionerParams ) );

    AMP::shared_ptr<AMP::Database> mechanicsPreconditioner_db =
        columnPreconditioner_db->getDatabase( "mechanicsPreconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> mechanicsPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( mechanicsPreconditioner_db ) );
    mechanicsPreconditionerParams->d_pOperator = linearMechanicsOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearMechanicsPreconditioner(
        new AMP::Solver::TrilinosMLSolver( mechanicsPreconditionerParams ) );

    AMP::shared_ptr<AMP::Database> thermalPreconditioner_db =
        columnPreconditioner_db->getDatabase( "thermalPreconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( thermalPreconditioner_db ) );
    thermalPreconditionerParams->d_pOperator = linearThermalOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalPreconditioner(
        new AMP::Solver::TrilinosMLSolver( thermalPreconditionerParams ) );

    columnPreconditioner->append( linearMechanicsPreconditioner );
    columnPreconditioner->append( linearThermalPreconditioner );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // register the preconditioner with the Jacobian free Krylov solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver =
        nonlinearSolver->getKrylovSolver();

    linearSolver->setPreconditioner( columnPreconditioner );

    nonlinearThermoMechanicsOperator->residual( rhsVec, solVec, resVec );
    double initialResidualNorm = resVec->L2Norm();

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;

    nonlinearSolver->setZeroInitialGuess( false );

    nonlinearSolver->solve( rhsVec, solVec );

    nonlinearThermoMechanicsOperator->residual( rhsVec, solVec, resVec );

    double finalResidualNorm = resVec->L2Norm();

    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    AMP::LinearAlgebra::Vector::shared_ptr mechUvec =
        displacementVec->select( AMP::LinearAlgebra::VS_Stride( 0, 3 ), "U" );
    AMP::LinearAlgebra::Vector::shared_ptr mechVvec =
        displacementVec->select( AMP::LinearAlgebra::VS_Stride( 1, 3 ), "V" );
    AMP::LinearAlgebra::Vector::shared_ptr mechWvec =
        displacementVec->select( AMP::LinearAlgebra::VS_Stride( 2, 3 ), "W" );

    double finalMaxU = mechUvec->maxNorm();
    double finalMaxV = mechVvec->maxNorm();
    double finalMaxW = mechWvec->maxNorm();

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


int main( int argc, char *argv[] )
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
