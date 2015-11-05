#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "operators/libmesh/VolumeIntegralOperator.h"
#include "operators/NeutronicsRhs.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "materials/Material.h"

#include "utils/Writer.h"

#include "ampmesh/Mesh.h"
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "operators/flow/NavierStokesLSWFLinearFEOperator.h"
#include "operators/flow/NavierStokesLSWFFEOperator.h"

#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/ColumnOperator.h"
#include "operators/OperatorBuilder.h"

#include "solvers/ColumnSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscSNESSolverParameters.h"
#include "solvers/petsc/PetscSNESSolver.h"

#include "solvers/trilinos/TrilinosMLSolver.h"


void myTest(AMP::UnitTest *ut, std::string exeName)
{
    std::string input_file = "input_" + exeName;
    //std::string log_file = "output_" + exeName;

    //  AMP::PIO::logOnlyNodeZero(log_file);
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

    size_t N_error0 = ut->NumFailLocal();;

    AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
    AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
    input_db->printClassData(AMP::plog);

    //--------------------------------------------------
    //   Create the Mesh.
    //--------------------------------------------------
    AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
    AMP::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
    mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh(mgrParams);
    //--------------------------------------------------

    //--------------------------------------------------
    // Create a DOF manager for a nodal vector 
    //--------------------------------------------------
    int DOFsPerNode = 10;
    int DOFsPerElement = 8;
    int nodalGhostWidth = 1;
    int gaussPointGhostWidth = 1;
    bool split = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap = 
        AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth,      DOFsPerNode,    split);
    AMP::Discretization::DOFManager::shared_ptr gaussPointDofMap = 
        AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Volume, gaussPointGhostWidth, DOFsPerElement, split);

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear flow 
    AMP_INSIST( input_db->keyExists("NonlinearFlowOperator"), "key missing!" );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> flowTransportModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearFlowOperator = 
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
			"NonlinearFlowOperator", input_db, flowTransportModel));

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the input variable
    AMP::shared_ptr<AMP::Operator::NavierStokesLSWFFEOperator> flowVolumeOperator =
	    AMP::dynamic_pointer_cast<AMP::Operator::NavierStokesLSWFFEOperator>(nonlinearFlowOperator->getVolumeOperator());

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> flowVariable = flowVolumeOperator->getOutputVariable();

    // create solution, rhs, and residual vectors
    AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( nodalDofMap, flowVariable);
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, flowVariable);
    AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( nodalDofMap, flowVariable);

    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for flow
    AMP_INSIST( input_db->keyExists("LinearFlowOperator"), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearFlowOperator = 
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator( meshAdapter,
		"LinearFlowOperator", input_db, flowTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    //Initial guess

    solVec->setToScalar(0.);
    double initialGuessNorm  = solVec->L2Norm();
    std::cout << "initial guess norm = " << initialGuessNorm <<"\n";

    nonlinearFlowOperator->modifyInitialSolutionVector(solVec);

    initialGuessNorm  = solVec->L2Norm();
    std::cout << "initial guess norm  after apply = " << initialGuessNorm <<"\n";

    rhsVec->zero();

    nonlinearFlowOperator->modifyRHSvector(rhsVec);

    double initialRhsNorm  = rhsVec->L2Norm();
    std::cout << "rhs norm  after modifyRHSvector = " << initialRhsNorm <<"\n";
    double expectedVal = 0.;
    if( !AMP::Utilities::approx_equal( expectedVal, initialRhsNorm, 1e-5) )
        ut->failure("the rhs norm after modifyRHSvector has changed.");

    // Get the solver databases
    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
    AMP::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

    // Create the preconditioner
    AMP::shared_ptr<AMP::Database> flowPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> flowPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters(flowPreconditioner_db));
    flowPreconditionerParams->d_pOperator = linearFlowOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFlowPreconditioner(new AMP::Solver::TrilinosMLSolver(flowPreconditionerParams));

    // initialize the linear solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
        AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
    linearSolverParams->d_pOperator = linearFlowOperator;
    linearSolverParams->d_comm = globalComm;
    linearSolverParams->d_pPreconditioner = linearFlowPreconditioner;
//    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));

    // Crete the solvers
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
        new AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));
    nonlinearSolverParams->d_comm = globalComm;
    nonlinearSolverParams->d_pOperator = nonlinearFlowOperator;
//    nonlinearSolverParams->d_pKrylovSolver = linearSolver;
    nonlinearSolverParams->d_pInitialGuess = solVec;
    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();
    linearSolver->setPreconditioner(linearFlowPreconditioner);

    nonlinearFlowOperator->residual(rhsVec, solVec, resVec);
    double initialResidualNorm  = resVec->L2Norm();

    AMP::pout<<"Initial Residual Norm: "<<initialResidualNorm<<std::endl;
    expectedVal = 3625.84;
    if( !AMP::Utilities::approx_equal( expectedVal, initialResidualNorm, 1e-5) ) {
        ut->failure("the Initial Residual Norm has changed."); }

    nonlinearSolver->setZeroInitialGuess(false);

    nonlinearSolver->solve(rhsVec, solVec);

    solVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    resVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

    nonlinearFlowOperator->residual(rhsVec, solVec, resVec);

    double finalResidualNorm  = resVec->L2Norm();
    double finalSolutionNorm  = solVec->L2Norm();
    double finalRhsNorm  = rhsVec->L2Norm();

    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;
    std::cout << "Final Solution Norm: " << finalSolutionNorm << std::endl;
    std::cout << "Final Rhs Norm: "      << finalRhsNorm      << std::endl;

    if( fabs(finalResidualNorm) > 1e-9 )
        ut->failure("the Final Residual is larger than the tolerance");
    if( !AMP::Utilities::approx_equal( 45431.3, solVec->L2Norm(), 1e-5) )
        ut->failure("the Final Solution Norm has changed.");
    if( !AMP::Utilities::approx_equal( initialRhsNorm, finalRhsNorm, 1e-12) )
        ut->failure("the Final Rhs Norm has changed.");

    #ifdef USE_EXT_SILO
        AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
        siloWriter->registerMesh( meshAdapter );
        siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
        siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::Vertex, "Residual" );
        siloWriter->writeFile( input_file , 0 );
    #endif

    if ( N_error0 == ut->NumFailLocal() )
        ut->passes(exeName);
    else
        ut->failure(exeName);

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back("testPetscSNESSolver-IncompressibleFlow-1");

    for(unsigned int i = 0; i < exeNames.size(); i++) {
        try {
            myTest(&ut, exeNames[i]);
        } catch (std::exception &err) {
            std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
            ut.failure("ERROR: While testing");
        } catch( ... ) {
            std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
            ut.failure("ERROR: While testing");
        }
    }
   
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   

