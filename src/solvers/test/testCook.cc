#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include <cassert>
#include <fstream>

#include <sys/stat.h>

/* Boost files */
#include "boost/shared_ptr.hpp"

/* libMesh files */
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/elem.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/boundary_info.h"
#include "libmesh/fe_base.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"

/* AMP files */
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/initializeLibMesh.h"
#include "ampmesh/libmesh/libMesh.h"
#include "materials/Material.h"

#include "operators/LinearBVPOperator.h"
#include "operators/ElementPhysicsModel.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/DirichletVectorCorrection.h"


#include "vectors/Vector.h"
#include "vectors/VectorSelector.h"
#include "utils/Writer.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"

#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

#include "utils/ReadTestMesh.h"

extern "C"{
#include "petsc.h"
}


void linearElasticTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName + ".txt";

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  boost::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(new AMP::Mesh::initializeLibMesh(globalComm));

  const unsigned int mesh_dim = 3;
  boost::shared_ptr< ::Mesh > mesh(new ::Mesh(mesh_dim));

  std::string mesh_file = input_db->getString("mesh_file");
  if(globalComm.getRank() == 0) {
    AMP::readTestMesh(mesh_file, mesh);
  }//end if root processor

  MeshCommunication().broadcast(*(mesh.get()));
  mesh->prepare_for_use(false);
  AMP::Mesh::Mesh::shared_ptr meshAdapter ( new AMP::Mesh::libMesh(mesh,"cook") );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														 "MechanicsBVPOperator",
														 input_db,
														 elementPhysicsModel));

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
															 "Load_Boundary",
															 input_db,
															 dummyModel));
  //This has an in-place apply. So, it has an empty input variable and
  //the output variable is the same as what it is operating on. 
  dirichletVecOp->setVariable(bvpOperator->getOutputVariable());

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  AMP::Discretization::DOFManager::shared_ptr DOF_vector = AMP::Discretization::simpleDOFManager::create(meshAdapter,AMP::Mesh::Vertex,1,3,true);
  AMP::LinearAlgebra::Vector::shared_ptr mechSolVec = AMP::LinearAlgebra::createVector( DOF_vector, bvpOperator->getOutputVariable(), true );
  AMP::LinearAlgebra::Vector::shared_ptr mechRhsVec = AMP::LinearAlgebra::createVector( DOF_vector, bvpOperator->getOutputVariable(), true );
  AMP::LinearAlgebra::Vector::shared_ptr mechResVec = AMP::LinearAlgebra::createVector( DOF_vector, bvpOperator->getOutputVariable(), true );

  mechSolVec->setToScalar(0.5);
  mechRhsVec->setToScalar(0.0);
  mechResVec->setToScalar(0.0);

  dirichletVecOp->apply(nullVec, nullVec, mechRhsVec, 1.0, 0.0);

  double rhsNorm = mechRhsVec->L2Norm();

  AMP::pout<<"RHS Norm: "<<rhsNorm<<std::endl;

  double initSolNorm = mechSolVec->L2Norm();

  AMP::pout<<"Initial Solution Norm: "<<initSolNorm<<std::endl;

  bvpOperator->apply(mechRhsVec, mechSolVec, mechResVec, 1.0, -1.0);

  double initResidualNorm = mechResVec->L2Norm();

  AMP::pout<<"Initial Residual Norm: "<<initResidualNorm<<std::endl;

  boost::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase("LinearSolver"); 

  // ---- first initialize the preconditioner
  boost::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase("Preconditioner"); 
  boost::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(new AMP::Solver::TrilinosMLSolverParameters(pcSolver_db));
  pcSolverParams->d_pOperator = bvpOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(new AMP::Solver::TrilinosMLSolver(pcSolverParams));

  // initialize the linear solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = bvpOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = pcSolver;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));

  linearSolver->setZeroInitialGuess(false);

  linearSolver->solve(mechRhsVec, mechSolVec);

  double finalSolNorm = mechSolVec->L2Norm();

  AMP::pout<<"Final Solution Norm: "<<finalSolNorm<<std::endl;

  AMP::LinearAlgebra::Vector::shared_ptr mechUvec = mechSolVec->select( AMP::LinearAlgebra::VS_Stride(0,3), "U" );
  AMP::LinearAlgebra::Vector::shared_ptr mechVvec = mechSolVec->select( AMP::LinearAlgebra::VS_Stride(1,3), "V" );
  AMP::LinearAlgebra::Vector::shared_ptr mechWvec = mechSolVec->select( AMP::LinearAlgebra::VS_Stride(2,3), "W" );

  double finalMaxU = mechUvec->maxNorm();
  double finalMaxV = mechVvec->maxNorm();
  double finalMaxW = mechWvec->maxNorm();

  AMP::pout<<"Maximum U displacement: "<<finalMaxU<<std::endl;
  AMP::pout<<"Maximum V displacement: "<<finalMaxV<<std::endl;
  AMP::pout<<"Maximum W displacement: "<<finalMaxW<<std::endl;

  bvpOperator->apply(mechRhsVec, mechSolVec, mechResVec, 1.0, -1.0);

  double finalResidualNorm = mechResVec->L2Norm();

  AMP::pout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

  if(finalResidualNorm > (1e-10*initResidualNorm)) {
    ut->failure(exeName);
  } else {
    ut->passes(exeName);
  }

}

int main(int argc, char *argv[]) {

    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;

    if(argc == 1) {
        exeNames.push_back("testCook-normal-mesh0");
        exeNames.push_back("testCook-reduced-mesh0");

        exeNames.push_back("testCook-normal-mesh1");
        exeNames.push_back("testCook-reduced-mesh1");

        exeNames.push_back("testCook-normal-mesh2");
        exeNames.push_back("testCook-reduced-mesh2");
    } else {
        for(int i = 1; i < argc; i+= 2) {
          char inpName[100];
          sprintf(inpName, "testCook-%s-mesh%d", argv[i], atoi(argv[i+1]));
          exeNames.push_back(inpName);
        }//end for i
    }

    for(size_t i = 0; i < exeNames.size(); i++) {
        try {
            linearElasticTest(&ut, exeNames[i]);
        } catch (std::exception &err) {
            AMP::pout << "ERROR: " << err.what() << std::endl;
            ut.failure("ERROR: While testing");
        } catch( ... ) {
            AMP::pout << "ERROR: " << "An unknown exception was thrown." << std::endl;
            ut.failure("ERROR: While testing");
        }
    } //end for i

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}   


