
#include <iostream>
#include <string>

#include "materials/Material.h"

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/ReadTestMesh.h"
#include "utils/WriteSolutionToFile.h"

#include "ampmesh/libmesh/initializeLibMesh.h"
#include "ampmesh/libmesh/libMesh.h"
#include "utils/ReadTestMesh.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "mesh_communication.h"

#include "vectors/trilinos/EpetraVector.h"

#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

#include "ml_include.h"

int myMatVec(ML_Operator *data, int in_length, double in[], int out_length, double out[]) {

  AMP::Operator::LinearOperator * op = reinterpret_cast<AMP::Operator::LinearOperator *>(ML_Get_MyMatvecData(data));
  AMP::LinearAlgebra::Matrix::shared_ptr mat = op->getMatrix();

  AMP::LinearAlgebra::Vector::shared_ptr inVec = mat->getRightVector();
  AMP::LinearAlgebra::Vector::shared_ptr outVec = mat->getLeftVector();

  inVec->putRawData(in);

  mat->mult(inVec, outVec);

  double* outPtr = outVec->getRawDataBlock<double>();
  for(int i = 0; i < out_length; i++) {
    out[i] = outPtr[i];
  }

  return 0;
}

int myGetRow(ML_Operator *data, int N_requested_rows, int requested_rows[],
    int allocated_space, int columns[], double values[], int row_lengths[] ) {

  AMP::Operator::LinearOperator * op = reinterpret_cast<AMP::Operator::LinearOperator *>(ML_Get_MyGetrowData(data));
  AMP::LinearAlgebra::Matrix::shared_ptr mat = op->getMatrix();

  int spaceRequired = 0;
  int cnt = 0;
  for(int i = 0; i < N_requested_rows; i++) {
    int row = requested_rows[i];
    std::vector<unsigned int> cols;
    std::vector<double> vals;

    mat->getRowByGlobalID(row, cols, vals);
    spaceRequired += cols.size();

    if(allocated_space >= spaceRequired) {
      for(size_t j = 0; j < cols.size(); j++) {
        columns[cnt] = cols[j];
        values[cnt] = vals[j];
        cnt++;
      }
      row_lengths[i] = cols.size();
    } else {
      return 0;
    }
  }

  return 1;
}

void myTest(AMP::UnitTest *ut, std::string exeName, int type) {
  std::string input_file = "input_" + exeName;
  char log_file[200];
  sprintf(log_file, "output_%s_%d", exeName.c_str(), type);

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);
  std::string mesh_file = input_db->getString("mesh_file");

  boost::shared_ptr<AMP::InputDatabase> mesh_file_db(new AMP::InputDatabase("mesh_file_db"));
  AMP::InputManager::getManager()->parseInputFile(mesh_file, mesh_file_db);

  boost::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(new AMP::Mesh::initializeLibMesh(globalComm));

  const unsigned int mesh_dim = 3;
  boost::shared_ptr< ::Mesh > fusedMesh(new ::Mesh(mesh_dim));

  AMP::readTestMesh(mesh_file, fusedMesh);

  MeshCommunication().broadcast(*(fusedMesh.get()));

  fusedMesh->prepare_for_use(false);

  AMP::Mesh::Mesh::shared_ptr fusedMeshAdapter( new AMP::Mesh::libMesh( fusedMesh, "mesh" ) );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> fusedElementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> fusedOperator = boost::dynamic_pointer_cast<
  AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(fusedMeshAdapter,
										   "BVPOperator",
										   input_db,
										   fusedElementPhysicsModel));

  AMP::LinearAlgebra::Variable::shared_ptr fusedVar = fusedOperator->getOutputVariable();

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOperator =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
									  AMP::Operator::OperatorBuilder::createOperator(fusedMeshAdapter, "LoadOperator", input_db, dummyModel));
  loadOperator->setVariable(fusedVar);

  AMP::Discretization::DOFManager::shared_ptr NodalVectorDOF = 
     AMP::Discretization::simpleDOFManager::create(fusedMeshAdapter,AMP::Mesh::Vertex,1,3);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr fusedSolVec = AMP::LinearAlgebra::createVector(NodalVectorDOF,fusedVar);
  AMP::LinearAlgebra::Vector::shared_ptr fusedRhsVec = fusedSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr fusedResVec = fusedSolVec->cloneVector();

  fusedRhsVec->zero();
  loadOperator->apply(nullVec, nullVec, fusedRhsVec, 1.0, 0.0);

  boost::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase("MLoptions"); 

  std::cout<<std::endl;

  size_t localSize = fusedSolVec->getLocalSize();

  //Matrix-based
  if(type == 0) {
    ML_set_random_seed(123456);
    std::cout<<"Matrix-Based ML: "<<std::endl;
    fusedSolVec->zero();

    boost::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> mlSolverParams(new AMP::Solver::TrilinosMLSolverParameters(mlSolver_db));
    boost::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(new AMP::Solver::TrilinosMLSolver(mlSolverParams));

    AMP::LinearAlgebra::Matrix::shared_ptr mat = fusedOperator->getMatrix();

    AMP::LinearAlgebra::Matrix::shared_ptr matCopy = mat->cloneMatrix();
    matCopy->zero();
    matCopy->axpy(1.0, mat);

    mat->zero();
    mlSolver->registerOperator(fusedOperator);

    mat->axpy(1.0, matCopy);

    mlSolver->solve(fusedRhsVec, fusedSolVec);
    std::cout<<std::endl;
  }

  //Matrix-Free-1
  if(type == 1) {
    ML_set_random_seed(123456);
    std::cout<<"Matrix-Free ML Type-1: "<<std::endl;
    fusedSolVec->zero();

    ML_Comm *comm;
    ML_Comm_Create(&comm);
    ML_Comm_Set_UsrComm(comm, globalComm.getCommunicator() );

    ML_Operator *ml_op = ML_Operator_Create(comm);
    ML_Operator_Set_ApplyFuncData(ml_op, localSize, localSize, fusedOperator.get(), localSize, myMatVec, 0);
    ML_Operator_Set_Getrow(ml_op, localSize, myGetRow);

    Teuchos::ParameterList paramsList;
    ML_Epetra::SetDefaults("SA", paramsList);
    paramsList.set("ML output", mlSolver_db->getInteger("print_info_level") );
    paramsList.set("PDE equations", mlSolver_db->getInteger("PDE_equations") );
    paramsList.set("cycle applications", mlSolver_db->getInteger("max_iterations") );
    paramsList.set("max levels", mlSolver_db->getInteger("max_levels") );

    boost::shared_ptr<ML_Epetra::MultiLevelPreconditioner> mlSolver( new 
        ML_Epetra::MultiLevelPreconditioner(ml_op, paramsList));

    const ML_Aggregate* agg_obj = mlSolver->GetML_Aggregate();
    ML_Aggregate_Print(const_cast<ML_Aggregate*>(agg_obj));

    Epetra_Vector &fVec = (AMP::LinearAlgebra::EpetraVector::view ( fusedRhsVec ))->castTo<AMP::LinearAlgebra::EpetraVector>().getEpetra_Vector();
    Epetra_Vector &uVec = (AMP::LinearAlgebra::EpetraVector::view ( fusedSolVec ))->castTo<AMP::LinearAlgebra::EpetraVector>().getEpetra_Vector();

    fusedOperator->apply(fusedRhsVec, fusedSolVec, fusedResVec, 1.0, -1.0);
    std::cout << "MatFree-1: L2 norm of residual before solve " <<std::setprecision(15)<< fusedResVec->L2Norm() << std::endl;

    mlSolver->ApplyInverse(fVec, uVec);

    if ( fusedSolVec->isA<AMP::LinearAlgebra::DataChangeFirer>() )
    {
      fusedSolVec->castTo<AMP::LinearAlgebra::DataChangeFirer>().fireDataChange();
    }

    double solution_norm = fusedSolVec->L2Norm();
    std::cout << "MatFree-1:  solution norm: " <<std::setprecision(15)<< solution_norm << std::endl;

    fusedOperator->apply(fusedRhsVec, fusedSolVec, fusedResVec, 1.0, -1.0);
    std::cout << "MatFree-1: L2 norm of residual after solve " <<std::setprecision(15)<< fusedResVec->L2Norm() << std::endl;    

    ML_Operator_Destroy(&ml_op);

    ML_Comm_Destroy(&comm);
  }

  //Matrix-Free-2
  if(type == 2) {
    ML_set_random_seed(123456);
    std::cout<<"Matrix-Free ML Type-2: "<<std::endl;
    fusedSolVec->zero();

    int numGrids = mlSolver_db->getInteger("max_levels");
    int numPDEs = mlSolver_db->getInteger("PDE_equations");

    ML *ml_object;
    ML_Create (&ml_object, numGrids);

    ML_Init_Amatrix(ml_object, 0, localSize, localSize, fusedOperator.get());
    ML_Set_Amatrix_Getrow(ml_object, 0, &myGetRow, NULL, localSize);
    ML_Set_Amatrix_Matvec(ml_object, 0, &myMatVec);
    ML_Set_MaxIterations(ml_object, 1 + mlSolver_db->getInteger("max_iterations") );
    ML_Set_ResidualOutputFrequency(ml_object, 1);
    ML_Set_PrintLevel( mlSolver_db->getInteger("print_info_level") );
    ML_Set_OutputLevel(ml_object,  mlSolver_db->getInteger("print_info_level") );

    ML_Aggregate *agg_object;
    ML_Aggregate_Create(&agg_object);
    agg_object->num_PDE_eqns = numPDEs;
    agg_object->nullspace_dim = numPDEs;
    ML_Aggregate_Set_MaxCoarseSize(agg_object, 128);
    ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(agg_object);

    int nlevels = ML_Gen_MGHierarchy_UsingAggregation(ml_object, 0, ML_INCREASING, agg_object);
    std::cout<<"Number of actual levels : "<< nlevels <<std::endl;

    for(int lev = 0; lev < (nlevels - 1); lev++) {
      ML_Gen_Smoother_SymGaussSeidel(ml_object, lev, ML_BOTH, 2, 1.0);
    }
    ML_Gen_Smoother_Amesos(ml_object, (nlevels - 1), ML_AMESOS_KLU, -1, 0.0);

    ML_Gen_Solver(ml_object, ML_MGV, 0, (nlevels-1));

    double * solArr = fusedSolVec->getRawDataBlock<double>();
    double * rhsArr = fusedRhsVec->getRawDataBlock<double>();

    fusedOperator->apply(fusedRhsVec, fusedSolVec, fusedResVec, 1.0, -1.0);
    std::cout << "MatFree-2: L2 norm of residual before solve " <<std::setprecision(15)<< fusedResVec->L2Norm() << std::endl;

    ML_Iterate(ml_object, solArr, rhsArr);

    if ( fusedSolVec->isA<AMP::LinearAlgebra::DataChangeFirer>() )
    {
      fusedSolVec->castTo<AMP::LinearAlgebra::DataChangeFirer>().fireDataChange();
    }

    double solution_norm = fusedSolVec->L2Norm();
    std::cout << "MatFree-2:  solution norm: " <<std::setprecision(15)<< solution_norm << std::endl;

    fusedOperator->apply(fusedRhsVec, fusedSolVec, fusedResVec, 1.0, -1.0);
    std::cout << "MatFree-2: L2 norm of residual after solve " <<std::setprecision(15)<< fusedResVec->L2Norm() << std::endl;    

    ML_Aggregate_Destroy(&agg_object);
    ML_Destroy(&ml_object);
    std::cout<<std::endl;
  }

  char outFile[200];
  sprintf(outFile, "%s-%d", exeName.c_str(), type);
  printSolution(fusedMeshAdapter, fusedSolVec, outFile);

  ut->passes(exeName);
}

void loopMyTest(AMP::UnitTest *ut, std::string exeName) {
  for(int type = 0; type < 3; type++) {
    myTest(ut, exeName, type);
  }
}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::vector<std::string> exeNames;

  if(argc == 1) {
    exeNames.push_back("testMatrixFreeML-1");
  } else {
    for(int i = 1; i < argc; i++) {
      char inpName[100];
      sprintf(inpName, "testMatrixFreeML-%s", argv[i]);
      exeNames.push_back(inpName);
    }//end for i
  }

  for(size_t i = 0; i < exeNames.size(); i++) {
    try {
      loopMyTest(&ut, exeNames[i]);
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


