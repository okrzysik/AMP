#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

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

#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/ColumnOperator.h"
#include "operators/OperatorBuilder.h"

#include "../ColumnSolver.h"
#include "../PetscKrylovSolverParameters.h"
#include "../PetscKrylovSolver.h"
#include "../TrilinosMLSolver.h"


void myTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

//--------------------------------------------------
//   Create the Mesh.
//--------------------------------------------------
    AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
    boost::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
    boost::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
    mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
    boost::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh(mgrParams);
//--------------------------------------------------

//--------------------------------------------------
// Create a DOF manager for a nodal vector 
//--------------------------------------------------
  int DOFsPerNode = 1;
  int nodalGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr nodalDofMap      = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth,      DOFsPerNode,    split);
//--------------------------------------------------

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> FickMaterialModel;
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for Fick
  AMP_INSIST( input_db->keyExists("testLinearFickOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearFickOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
																	 AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																							"testLinearFickOperator",
																							input_db,
																							FickMaterialModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for thermal
  AMP_INSIST( input_db->keyExists("testLinearThermalOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
																	    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																							   "testLinearThermalOperator",
																							   input_db,
																							   thermalTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a column operator object for linear Thermal-Fick
  boost::shared_ptr<AMP::Operator::OperatorParameters> params;
  boost::shared_ptr<AMP::Operator::ColumnOperator> linearThermalFickOperator(new AMP::Operator::ColumnOperator(params));
  linearThermalFickOperator->append(linearThermalOperator);
  linearThermalFickOperator->append(linearFickOperator);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the input multi-variable
  boost::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> thermalVolumeOperator = boost::dynamic_pointer_cast<
  AMP::Operator::DiffusionLinearFEOperator>(linearThermalOperator->getVolumeOperator());
  boost::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> FickVolumeOperator = boost::dynamic_pointer_cast<
    AMP::Operator::DiffusionLinearFEOperator>(linearFickOperator->getVolumeOperator());

  boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputVariable(new AMP::LinearAlgebra::MultiVariable("inputVariable"));
  inputVariable->add(thermalVolumeOperator->getInputVariable());
  inputVariable->add(FickVolumeOperator->getInputVariable());

  // initialize the output multi-variable
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable = linearThermalFickOperator->getOutputVariable();

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( nodalDofMap, inputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, outputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( nodalDofMap, outputVariable );

  // create the following shared pointers for ease of use
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //Initial-Guess for thermal
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyThermalModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletThermalInVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
															 "ThermalInitialGuess",
															 input_db,
															 dummyThermalModel));
  dirichletThermalInVecOp->setVariable(thermalVolumeOperator->getInputVariable());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //Initial-Guess for Fick
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyFickModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletFickInVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
															 "FickInitialGuess",
															 input_db,
															 dummyFickModel));
  dirichletFickInVecOp->setVariable(FickVolumeOperator->getInputVariable());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //Random initial guess
  solVec->setRandomValues();

  //Initial guess for thermal must satisfy the thermal Dirichlet boundary conditions
  dirichletThermalInVecOp->apply(nullVec, nullVec, solVec, 1.0, 0.0);

  //Initial guess for Fick must satisfy the Fick Dirichlet boundary conditions
  dirichletFickInVecOp->apply(nullVec, nullVec, solVec, 1.0, 0.0);

  rhsVec->setToScalar(0.0);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // ---- initialize the solvers
  boost::shared_ptr<AMP::Database> columnSolver_db = input_db->getDatabase("ColumnSolver");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> columnSolverParams(new AMP::Solver::SolverStrategyParameters(columnSolver_db));
  boost::shared_ptr<AMP::Solver::ColumnSolver> columnSolver(new AMP::Solver::ColumnSolver(columnSolverParams));
    
  boost::shared_ptr<AMP::Database> thermalSolver_db = columnSolver_db->getDatabase("ThermalSolver"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalSolverParams(new AMP::Solver::SolverStrategyParameters(thermalSolver_db));
  thermalSolverParams->d_pOperator = linearThermalOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalSolver(new AMP::Solver::TrilinosMLSolver(thermalSolverParams));

  boost::shared_ptr<AMP::Database> FickSolver_db = columnSolver_db->getDatabase("FickSolver"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> FickSolverParams(new AMP::Solver::SolverStrategyParameters(FickSolver_db));
  FickSolverParams->d_pOperator = linearFickOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFickSolver(new AMP::Solver::TrilinosMLSolver(FickSolverParams));

  columnSolver->append(linearThermalSolver);
  columnSolver->append(linearFickSolver);
  
  double initialResidualNorm  = resVec->L2Norm();

  AMP::pout<<"Initial Residual Norm: "<<initialResidualNorm<<std::endl;

  columnSolver->solve(rhsVec, solVec);

  linearThermalFickOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);

  double finalResidualNorm  = resVec->L2Norm();

  std::cout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

  if(finalResidualNorm>1.0e-08) {
    ut->failure("ColumnSolver unsuccessfully solves two linear operators");
  } else {
    ut->passes("ColumnSolver successfully solves two linear operators");
  }
  ut->passes(exeName);

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back("testColumnSolver-LinearThermalAndFick-1");

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

