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
#include "utils/PIO.h"
#include "materials/Material.h"


#include "utils/Writer.h"

#include "ampmesh/Mesh.h"
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"


#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"

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
    std::string log_file = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero(log_file);
    AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);

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
    int DOFsPerNode = 1;
    int DOFsPerElement = 8;
    int nodalGhostWidth = 1;
    int gaussPointGhostWidth = 1;
    bool split = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap      = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth,      DOFsPerNode,    split);
    AMP::Discretization::DOFManager::shared_ptr gaussPointDofMap = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Volume, gaussPointGhostWidth, DOFsPerElement, split);
    //--------------------------------------------------

    AMP::pout<<"Constructing Nonlinear Thermal Operator..."<<std::endl;

    //-------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear thermal diffusion
    AMP_INSIST( input_db->keyExists("testNonlinearThermalOperator"), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator = AMP::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
										    "testNonlinearThermalOperator",
										    input_db,
										    thermalTransportModel));

    //-------------------------------------------------------------------------------------------//
    // initialize the input variable
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator =
    AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearThermalOperator->getVolumeOperator());

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> thermalVariable = thermalVolumeOperator->getOutputVariable();

    // create solution, rhs, and residual vectors
    AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );

    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    //-------------------------------------------------------------------------------------------//

    AMP::pout<<"Constructing Linear Thermal Operator..."<<std::endl;

    //-------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists("testLinearThermalOperator"), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
																        AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																						       "testLinearThermalOperator",
																						       input_db,
																						       thermalTransportModel));

    ////////////////////////////////////
    //  CREATE THE NEUTRONICS SOURCE  //
    ////////////////////////////////////
    AMP_INSIST(input_db->keyExists("NeutronicsOperator"), "Key ''NeutronicsOperator'' is missing!");
    AMP::shared_ptr<AMP::Database>  neutronicsOp_db = input_db->getDatabase("NeutronicsOperator");
    AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ));
    AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(new AMP::Operator::NeutronicsRhs( neutronicsParams ));

    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar = neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr   SpecificPowerVec = AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

    neutronicsOperator->apply(nullVec, SpecificPowerVec);

    /////////////////////////////////////////////////////
    //  Integrate Nuclear Rhs over Desnity * Volume //
    /////////////////////////////////////////////////////

    AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator = AMP::dynamic_pointer_cast<
    AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
										      "VolumeIntegralOperator",
										      input_db,
										      stransportModel));

    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr   PowerInWattsVec = AMP::LinearAlgebra::createVector( nodalDofMap, PowerInWattsVar );
    PowerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply(SpecificPowerVec, PowerInWattsVec);

    rhsVec->copyVector(PowerInWattsVec);

    AMP::pout << "RHS L2 norm before corrections = " << std::setprecision(10) << (rhsVec->L2Norm()) <<"\n";
    AMP::pout << "RHS max before corrections = " << std::setprecision(10)<< (rhsVec->max()) <<"\n";
    AMP::pout << "RHS min before corrections = " << std::setprecision(10)<< (rhsVec->min()) <<"\n";

    nonlinearThermalOperator->modifyRHSvector(rhsVec);

    AMP::pout << "RHS L2 norm after corrections = " << std::setprecision(10) << (rhsVec->L2Norm()) <<"\n";
    AMP::pout << "RHS max after corrections = " << std::setprecision(10) << (rhsVec->max()) <<"\n";
    AMP::pout << "RHS min after corrections = " << std::setprecision(10) << (rhsVec->min()) <<"\n";

    //---------------------------------------------------------------------------------------------//
    //Initial guess

    double initGuess = input_db->getDoubleWithDefault("InitialGuess", 400.0);
    solVec->setToScalar(initGuess);

    AMP::pout << "initial guess L2 norm before corrections = " << (solVec->L2Norm()) <<"\n";
    AMP::pout << "initial guess max before corrections = " << (solVec->max()) <<"\n";
    AMP::pout << "initial guess min before corrections = " << (solVec->min()) <<"\n";

    nonlinearThermalOperator->modifyInitialSolutionVector(solVec);

    AMP::pout << "initial guess L2 norm after corrections = " << (solVec->L2Norm()) <<"\n";
    AMP::pout << "initial guess max after corrections = " << (solVec->max()) <<"\n";
    AMP::pout << "initial guess min after corrections = " << (solVec->min()) <<"\n";

    //---------------------------------------------------------------------------------------------/

    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
    AMP::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

    //---------------------------------------------------------------------------------------------//
    // initialize the nonlinear solver
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
        AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm = globalComm;
    nonlinearSolverParams->d_pOperator = nonlinearThermalOperator;
    nonlinearSolverParams->d_pInitialGuess = solVec;

    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

    AMP::pout<<"Calling Get Jacobian Parameters..."<<std::endl;

    //-------------------------------------------------------------------------------------------//
    nonlinearThermalOperator->modifyInitialSolutionVector(solVec);
    linearThermalOperator->reset(nonlinearThermalOperator->getParameters("Jacobian", solVec));

    AMP::pout<<"Finished reseting the jacobian."<<std::endl;

    //------------------------------------------------------------------------------------------//
    AMP::shared_ptr<AMP::Database> thermalPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams(new AMP::Solver::SolverStrategyParameters(thermalPreconditioner_db));
    thermalPreconditionerParams->d_pOperator = linearThermalOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalPreconditioner(new AMP::Solver::TrilinosMLSolver(thermalPreconditionerParams));

    //--------------------------------------------------------------------------------------//
    // register the preconditioner with the Jacobian free Krylov solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();
    linearSolver->setPreconditioner(linearThermalPreconditioner);
    nonlinearThermalOperator->residual(rhsVec, solVec, resVec);

    double initialResidualNorm  = resVec->L2Norm();
    AMP::pout<<"Initial Residual Norm: "<< initialResidualNorm <<std::endl;

    double expectedVal =20.7018 ;
    if( !AMP::Utilities::approx_equal( expectedVal, initialResidualNorm, 1e-5) ) {
        ut->failure("the Initial Residual Norm has changed."); }

    nonlinearSolver->setZeroInitialGuess(false);
    nonlinearSolver->solve(rhsVec, solVec);

    std::cout<<"Final Solution Norm: "<<solVec->L2Norm()<<std::endl;
    expectedVal = 45612 ;
    if( !AMP::Utilities::approx_equal( expectedVal, solVec->L2Norm(), 1e-5) ) {
        ut->failure("the Final Solution Norm has changed."); }

    AMP::pout<<" Solution Max: "<< std::setprecision(10) << (solVec->max()) <<std::endl;
    AMP::pout<<" Solution Min: "<< std::setprecision(10) <<(solVec->min()) <<std::endl;
    AMP::pout<<" Solution L1 Norm: "<< std::setprecision(10) <<(solVec->L1Norm()) <<std::endl;
    AMP::pout<<" Solution L2 Norm: "<< std::setprecision(10) <<(solVec->L2Norm()) <<std::endl;

    nonlinearThermalOperator->residual(rhsVec, solVec, resVec);

    double finalResidualNorm  = resVec->L2Norm();
    double finalSolutionNorm  = solVec->L2Norm();
    AMP::pout<<"Final Residual Norm: "<< std::setprecision(10) <<finalResidualNorm <<std::endl;
    AMP::pout<<"Final Solution Norm: "<< std::setprecision(10) <<finalSolutionNorm <<std::endl;

    expectedVal = 4.561204386863e4;
    if( fabs(finalResidualNorm) > 1e-8 )
        ut->failure("the Final Residual is larger than the tolerance");
    if( !AMP::Utilities::approx_equal( expectedVal, finalSolutionNorm, 1e-7) ) {
        ut->failure("the Final Residual Norm has changed."); }

    #ifdef USE_EXT_SILO
        AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
        siloWriter->registerMesh( meshAdapter );
        siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
        siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::Vertex, "Residual" );
        siloWriter->writeFile( exeName , 0 );
    #endif

    ut->passes(exeName);
}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back("testPetscSNESSolver-NonlinearThermal-cylinder_MATPRO2");

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

