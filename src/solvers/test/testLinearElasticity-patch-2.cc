
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include <iostream>
#include <string>
#include <cassert>
#include <fstream>

/* libMesh files */
#include "mesh.h"
#include "mesh_generation.h"
#include "mesh_communication.h"
#include "elem.h"
#include "cell_hex8.h"
#include "boundary_info.h"

/* AMP files */
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/libMesh.h"

#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "vectors/Vector.h"
#include "utils/Writer.h"

#include "solvers/trilinos/TrilinosMLSolver.h"

#include "utils/ReadTestMesh.h"

void linearElasticTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);

  if (globalComm.getSize() == 1) {
    boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
    AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
    input_db->printClassData(AMP::plog);

    boost::shared_ptr<AMP::InputDatabase> mesh_file_db(new AMP::InputDatabase("mesh_file_db"));
    AMP::InputManager::getManager()->parseInputFile((input_db->getString("mesh_file")), mesh_file_db);

    const unsigned int mesh_dim = 3;
    boost::shared_ptr< ::Mesh > mesh(new ::Mesh(mesh_dim));

    AMP::readTestMesh(mesh_file_db, mesh);
    mesh->prepare_for_use(false);

    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::shared_ptr (
        new AMP::Mesh::libMesh (mesh, "TestMesh") );

    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
      boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
            "MechanicsBVPOperator", input_db, elementPhysicsModel));

    AMP::LinearAlgebra::Variable::shared_ptr var = bvpOperator->getOutputVariable();

    AMP::Discretization::DOFManager::shared_ptr dofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::Vertex, 1, 3, true); 

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Vector::shared_ptr mechSolVec = AMP::LinearAlgebra::createVector(dofMap, var, true);
    AMP::LinearAlgebra::Vector::shared_ptr mechResVec = mechSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr mechRhsVec = mechSolVec->cloneVector();

    mechSolVec->setToScalar(0.5);
    mechResVec->setToScalar(0.0);
    mechRhsVec->setToScalar(0.0);

    bvpOperator->modifyRHSvector(mechRhsVec);

    double rhsNorm = mechRhsVec->L2Norm();

    AMP::pout<<"RHS Norm: "<<rhsNorm<<std::endl;

    double initSolNorm = mechSolVec->L2Norm();

    AMP::pout<<"Initial Solution Norm: "<<initSolNorm<<std::endl;

    bvpOperator->apply(mechRhsVec, mechSolVec, mechResVec, 1.0, -1.0);

    double initResidualNorm = mechResVec->L2Norm();

    AMP::pout<<"Initial Residual Norm: "<<initResidualNorm<<std::endl;

    boost::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase("LinearSolver"); 

    boost::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams(new AMP::Solver::SolverStrategyParameters(mlSolver_db));

    mlSolverParams->d_pOperator = bvpOperator;

    // create the ML solver interface
    boost::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(new AMP::Solver::TrilinosMLSolver(mlSolverParams));

    mlSolver->setZeroInitialGuess(false);

    mlSolver->solve(mechRhsVec, mechSolVec);

    double finalSolNorm = mechSolVec->L2Norm();

    AMP::pout<<"Final Solution Norm: "<<finalSolNorm<<std::endl;

    std::string fname = exeName + "-stressAndStrain.txt";

    (boost::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(bvpOperator->
                                                                           getVolumeOperator()))->printStressAndStrain(mechSolVec, fname);

    bvpOperator->apply(mechRhsVec, mechSolVec, mechResVec, 1.0, -1.0);

    double finalResidualNorm = mechResVec->L2Norm();

    AMP::pout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

    if(finalResidualNorm > (1e-10*initResidualNorm)) {
      ut->failure(exeName);
    } else {
      ut->passes(exeName);
    }
  } else {
    AMP::pout << "WARNING: This is a single processor test!" << std::endl;
    ut->passes(exeName);
  }
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  boost::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(new AMP::Mesh::initializeLibMesh(AMP_COMM_WORLD));

  AMP::UnitTest ut;

  std::vector<std::string> exeNames;
  exeNames.push_back("testLinearElasticity-patch-2-normal");
  exeNames.push_back("testLinearElasticity-patch-2-reduced");

  for(size_t i = 0; i < exeNames.size(); i++) {
    try {
      linearElasticTest(&ut, exeNames[i]);
    } catch (std::exception &err) {
      AMP::pout << "ERROR: " << err.what() << std::endl;
      ut.failure("ERROR");
    } catch( ... ) {
      AMP::pout << "ERROR: " << "An unknown exception was thrown." << std::endl;
      ut.failure("ERROR");
    }//end for reduced
  } //end for i

  ut.report();
  int num_failed = ut.NumFailGlobal();

  libmeshInit.reset();
  AMP::AMPManager::shutdown();
  return num_failed;

}

