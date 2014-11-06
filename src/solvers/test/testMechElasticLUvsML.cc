
#include <iostream>

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "libmesh/mesh_communication.h"

#include "mpi.h"

#include "utils/shared_ptr.h"

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>


#include "operators/OperatorBuilder.h"
#include "operators/LinearBVPOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

#include "ampmesh/libmesh/initializeLibMesh.h"
#include "ampmesh/libmesh/libMesh.h"
#include "utils/ReadTestMesh.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"


void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  int numMeshes = input_db->getInteger("NumberOfMeshFiles");
  AMP::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(new AMP::Mesh::initializeLibMesh(globalComm));

  AMP::shared_ptr<AMP::Database> ml_db = input_db->getDatabase("ML_Solver"); 
  AMP::shared_ptr<AMP::Database> lu_db = input_db->getDatabase("LU_Solver"); 
  AMP::shared_ptr<AMP::Database> cg_db = input_db->getDatabase("CG_Solver"); 
  AMP::shared_ptr<AMP::Database> rich_db = input_db->getDatabase("Richardson_Solver"); 

  for(int meshId = 1; meshId <= numMeshes; meshId++) {
    std::cout<<"Working on mesh "<<meshId<<std::endl;

    char meshFileKey[200];
    sprintf(meshFileKey, "mesh%d", meshId);

    std::string meshFile = input_db->getString(meshFileKey);

    const unsigned int mesh_dim = 3;
    AMP::shared_ptr< ::Mesh > mesh(new ::Mesh(mesh_dim));

    if(globalComm.getRank() == 0) {
      AMP::readBinaryTestMesh(meshFile, mesh);
    }

    MeshCommunication().broadcast(*(mesh.get()));
    mesh->prepare_for_use(false);

    AMP::Mesh::Mesh::shared_ptr meshAdapter( new AMP::Mesh::libMesh( mesh, "mesh" ) );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator = AMP::dynamic_pointer_cast<
      AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
            "BVPOperator", input_db, elementPhysicsModel));

    AMP::LinearAlgebra::Variable::shared_ptr dispVar = bvpOperator->getOutputVariable();

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOperator = AMP::dynamic_pointer_cast<
      AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
            "LoadOperator", input_db, dummyModel));
    loadOperator->setVariable(dispVar);

    AMP::Discretization::DOFManager::shared_ptr NodalVectorDOF = 
        AMP::Discretization::simpleDOFManager::create(meshAdapter,AMP::Mesh::Vertex,1,3);

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector(NodalVectorDOF,dispVar);
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector(NodalVectorDOF,dispVar);
    AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector(NodalVectorDOF,dispVar);

    rhsVec->zero();
    loadOperator->apply(nullVec, nullVec, rhsVec, 1.0, 0.0);

    solVec->zero();
    resVec->zero();

    size_t numDofs = solVec->getGlobalSize();

    if(globalComm.getRank() == 0) {
      std::cout<<"Solving using LU"<<std::endl;
    }
    globalComm.barrier();
    double luStartTime = AMP::AMP_MPI::time();

    AMP::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> luParams(new AMP::Solver::TrilinosMLSolverParameters(lu_db));
    luParams->d_pOperator = bvpOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> luPC(new AMP::Solver::TrilinosMLSolver(luParams));

    AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> richParams(new AMP::Solver::PetscKrylovSolverParameters(rich_db));
    richParams->d_pOperator = bvpOperator;
    richParams->d_comm = globalComm;
    richParams->d_pPreconditioner = luPC;
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> richSolver(new AMP::Solver::PetscKrylovSolver(richParams));
    richSolver->setZeroInitialGuess(true);

    richSolver->solve(rhsVec, solVec);

    globalComm.barrier();
    double luEndTime = AMP::AMP_MPI::time();

    PetscInt richIters = 0;
    KSP richKsp = richSolver->getKrylovSolver();
    KSPGetIterationNumber(richKsp, &richIters);
    AMP_INSIST(richIters <= 1, "Should not need more than 1 LU-Richardson iteration.");

    KSPConvergedReason richReason;
    KSPGetConvergedReason(richKsp, &richReason);
    AMP_INSIST( ( (richReason == KSP_CONVERGED_RTOL) ||
          (richReason == KSP_CONVERGED_ATOL) ), "KSP did not converge properly." );

    richSolver.reset();
    luPC.reset();

    solVec->zero();
    resVec->zero();

    if(globalComm.getRank() == 0) {
      std::cout<<"Solving using ML"<<std::endl;
    }
    globalComm.barrier();
    double mlStartTime = AMP::AMP_MPI::time();

    AMP::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> mlParams(new AMP::Solver::TrilinosMLSolverParameters(ml_db));
    mlParams->d_pOperator = bvpOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> mlPC(new AMP::Solver::TrilinosMLSolver(mlParams));

    AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> cgParams(new AMP::Solver::PetscKrylovSolverParameters(cg_db));
    cgParams->d_pOperator = bvpOperator;
    cgParams->d_comm = globalComm;
    cgParams->d_pPreconditioner = mlPC;
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> cgSolver(new AMP::Solver::PetscKrylovSolver(cgParams));
    cgSolver->setZeroInitialGuess(true);

    cgSolver->solve(rhsVec, solVec);

    globalComm.barrier();
    double mlEndTime = AMP::AMP_MPI::time();

    PetscInt cgIters = 0;
    KSP cgKsp = cgSolver->getKrylovSolver();
    KSPGetIterationNumber(cgKsp, &cgIters);

    KSPConvergedReason cgReason;
    KSPGetConvergedReason(cgKsp, &cgReason);
    AMP_INSIST( ( (cgReason == KSP_CONVERGED_RTOL) ||
          (cgReason == KSP_CONVERGED_ATOL) ), "KSP did not converge properly." );

    cgSolver.reset();
    mlPC.reset();

    if(globalComm.getRank() == 0) {
      std::cout<<"Result: "<<numDofs<<" & "<<globalComm.getSize()<<
        " & "<<cgIters<<" & "<<(luEndTime - luStartTime)<<" & "
        <<(mlEndTime - mlStartTime)<<" \\\\ "<<std::endl; 
    }

  }//end for meshId

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testMechElasticLUvsML";

  try {
    myTest(&ut, exeName);
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  ut.report();

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;
}  


