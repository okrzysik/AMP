
#include <iostream>
#include <string>
#include <cassert>
#include <fstream>

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
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/PIO.h"

#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/libmesh/PressureBoundaryOperator.h"

#include "operators/mechanics/MechanicsLinearFEOperator.h"

#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "vectors/Vector.h"
#include "vectors/VectorSelector.h"
#include "utils/Writer.h"

#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

#include "utils/ReadTestMesh.h"

void linearElasticTest(AMP::UnitTest *ut, std::string exeName, 
    int exampleNum) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName + ".txt";

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

#ifdef USE_EXT_SILO
  // Create the silo writer and register the data
  AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
#endif

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "MechanicsBVPOperator", input_db, elementPhysicsModel));

  AMP::LinearAlgebra::Variable::shared_ptr dispVar = bvpOperator->getOutputVariable();

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "Load_Boundary", input_db, dummyModel));
  //This has an in-place apply. So, it has an empty input variable and
  //the output variable is the same as what it is operating on. 
  dirichletVecOp->setVariable(dispVar);

  //Pressure RHS
  boost::shared_ptr<AMP::Operator::PressureBoundaryOperator> pressureLoadVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::PressureBoundaryOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "Pressure_Boundary", input_db, dummyModel));
  //This has an in-place apply. So, it has an empty input variable and
  //the output variable is the same as what it is operating on. 

  AMP::Discretization::DOFManager::shared_ptr dofMap = AMP::Discretization::simpleDOFManager::create(
      meshAdapter, AMP::Mesh::Vertex, 1, 3, true); 

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr mechSolVec = AMP::LinearAlgebra::createVector(dofMap, dispVar, true);
  AMP::LinearAlgebra::Vector::shared_ptr mechRhsVec = mechSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr mechResVec = mechSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr mechPressureVec = mechSolVec->cloneVector();

  mechSolVec->setToScalar(0.0);
  mechRhsVec->setToScalar(0.0);
  mechResVec->setToScalar(0.0);
  mechPressureVec->setToScalar(0.0);

  dirichletVecOp->apply(nullVec, nullVec, mechRhsVec, 1.0, 0.0);

  double rhsNorm = mechRhsVec->L2Norm();
  AMP::pout<<"RHS Norm after Dirichlet Apply: "<<rhsNorm<<std::endl;

  double pressNorm = mechPressureVec->L2Norm();
  AMP::pout<<"Pressure Norm before Apply: "<<pressNorm<<std::endl;

  //Applying the pressure load
  pressureLoadVecOp->addRHScorrection(mechPressureVec);

  pressNorm = mechPressureVec->L2Norm();
  AMP::pout<<"Pressure Norm after Apply: "<<pressNorm<<std::endl;

  mechRhsVec->add(mechRhsVec, mechPressureVec);

  rhsNorm = mechRhsVec->L2Norm();
  AMP::pout<<"Total RHS Norm: "<<rhsNorm<<std::endl;

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

  std::string fname = exeName + "_StressAndStrain.txt";

  (boost::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(bvpOperator->
                                                                         getVolumeOperator()))->printStressAndStrain(mechSolVec, fname);

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

  double epsilon = 1.0e-13*(((bvpOperator->getMatrix())->extractDiagonal())->L1Norm());
  AMP::pout<<"epsilon = "<<epsilon<<std::endl;

#ifdef USE_EXT_SILO
  siloWriter->registerVector(mechSolVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
  char outFileName1[256];
  sprintf(outFileName1, "undeformedBeam_%d", exampleNum);
  siloWriter->writeFile(outFileName1, 1);
  meshAdapter->displaceMesh(mechSolVec);
  char outFileName2[256];
  sprintf(outFileName2, "deformedBeam_%d", exampleNum);
  siloWriter->writeFile(outFileName2, 1);
#endif

}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::vector<std::string> exeNames;

  if(argc == 1) {
    exeNames.push_back("testLinearMechanics-PressureBoundary-1");
    exeNames.push_back("testLinearMechanics-PressureBoundary-HaldenPellet");
  } else {
    for(int i = 1; i < argc; ++i) {
      char inpName[100];
      sprintf(inpName, "testLinearMechanics-PressureBoundary-%s", argv[i]);
      exeNames.push_back(inpName);
    }//end for i
  }

  for(size_t i = 0; i < exeNames.size(); ++i) {
    linearElasticTest(&ut, exeNames[i], i);
  }

  ut.report();

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;
}  


