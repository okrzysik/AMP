
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "externVars.h"

#include "ampmesh/libmesh/libMesh.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/SiloIO.h"
#include "ampmesh/dendro/DendroSearch.h"

#include "operators/OperatorBuilder.h"
#include "operators/LinearBVPOperator.h"
#include "operators/ColumnOperator.h"
#include "operators/PetscMatrixShellOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/contact/NodeToSegmentConstraintsOperator.h"

#include "solvers/ColumnSolver.h"
#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/contact/MPCSolver.h"

#include "utils/ReadTestMesh.h"

#include <fstream>
#include <boost/lexical_cast.hpp>

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName; 

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

#ifdef USE_SILO
  AMP::Mesh::SiloIO::shared_ptr siloWriter(new AMP::Mesh::SiloIO);
#endif

  int npes = globalComm.getSize();
  int rank = globalComm.getRank();
  std::fstream fout;
  std::string fileName = "debug_driver_" + boost::lexical_cast<std::string>(rank);
  fout.open(fileName.c_str(), std::fstream::out);

  // Load the input file
  globalComm.barrier();
  double inpReadBeginTime = MPI_Wtime();

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  globalComm.barrier();
  double inpReadEndTime = MPI_Wtime();
  if(!rank) {
    std::cout<<"Finished parsing the input file in "<<(inpReadEndTime - inpReadBeginTime)<<" seconds."<<std::endl;
  }

  // Load the meshes
  globalComm.barrier();
  double meshBeginTime = MPI_Wtime();

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(globalComm);
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);

  globalComm.barrier();
  double meshEndTime = MPI_Wtime();
  if(!rank) {
    std::cout<<"Finished reading the mesh in "<<(meshEndTime - meshBeginTime)<<" seconds."<<std::endl;
  }

  // Build the contact operator
  AMP_INSIST(input_db->keyExists("ContactOperator"), "Key ''ContactOperator'' is missing!");
  boost::shared_ptr<AMP::Database> contact_db = input_db->getDatabase("ContactOperator");
  std::vector<AMP::Mesh::MeshID> meshIDs = meshAdapter->getBaseMeshIDs();
  boost::shared_ptr<AMP::Operator::NodeToSegmentConstraintsOperatorParameters> 
      contactOperatorParams( new AMP::Operator::NodeToSegmentConstraintsOperatorParameters(contact_db) );
  AMP_INSIST(contact_db->keyExists("MasterMeshIndex"), "Key ''MasterMeshIndex'' is missing!");
  AMP_INSIST(contact_db->keyExists("SlaveMeshIndex"), "Key ''SlaveMeshIndex'' is missing!");
  AMP::Mesh::MeshID masterMeshID = meshIDs[contact_db->getInteger("MasterMeshIndex")];
  AMP::Mesh::MeshID slaveMeshID = meshIDs[contact_db->getInteger("SlaveMeshIndex")];
  contactOperatorParams->d_MasterMeshID = masterMeshID;
  contactOperatorParams->d_SlaveMeshID = slaveMeshID;
  AMP_INSIST(contact_db->keyExists("MasterBoundaryID"), "Key ''MasterBoundaryID'' is missing!");
  AMP_INSIST(contact_db->keyExists("SlaveBoundaryID"), "Key ''SlaveBoundaryID'' is missing!");
  contactOperatorParams->d_MasterBoundaryID = contact_db->getInteger("MasterBoundaryID");
  contactOperatorParams->d_SlaveBoundaryID = contact_db->getInteger("SlaveBoundaryID");
  
  int dofsPerNode = 3;
  int nodalGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr dofManager = AMP::Discretization::simpleDOFManager::create(meshAdapter,
      AMP::Mesh::Vertex, nodalGhostWidth, dofsPerNode, split);
  contactOperatorParams->d_DOFsPerNode = dofsPerNode;
  contactOperatorParams->d_DOFManager = dofManager;

  contactOperatorParams->d_GlobalComm = globalComm;
  contactOperatorParams->d_Mesh = meshAdapter;

  boost::shared_ptr<AMP::Operator::NodeToSegmentConstraintsOperator> 
      contactOperator( new AMP::Operator::NodeToSegmentConstraintsOperator(contactOperatorParams) );

  // TODO: RESET IN CONSTRUCTOR?
  contactOperator->reset(contactOperatorParams);

  // build a column operator and a column preconditioner
  boost::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator(new AMP::Operator::ColumnOperator(emptyParams));

  boost::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase("LinearSolver"); 
  boost::shared_ptr<AMP::Database> columnPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::ColumnSolverParameters> columnPreconditionerParams(new
      AMP::Solver::ColumnSolverParameters(columnPreconditioner_db));
  columnPreconditionerParams->d_pOperator = columnOperator;
  boost::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

  // build the master and slave operators
  AMP::Mesh::Mesh::shared_ptr masterMeshAdapter = meshAdapter->Subset(masterMeshID);
  if (masterMeshAdapter != NULL) {
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> masterOperator = boost::dynamic_pointer_cast<
        AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(masterMeshAdapter,
                                                                                         "MasterBVPOperator",
                                                                                         input_db,
                                                                                         masterElementPhysicsModel));
    columnOperator->append(masterOperator);

    boost::shared_ptr<AMP::Database> masterSolver_db = columnPreconditioner_db->getDatabase("MasterSolver"); 
    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> masterSolverParams(new
        AMP::Solver::PetscKrylovSolverParameters(masterSolver_db));
    masterSolverParams->d_pOperator = masterOperator;
    masterSolverParams->d_comm = masterMeshAdapter->getComm();
//    masterSolverParams->d_comm = globalComm;
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> masterSolver(new AMP::Solver::PetscKrylovSolver(masterSolverParams));
    columnPreconditioner->append(masterSolver);
  } // end if

  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> slaveLoadOperator;
  AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter = meshAdapter->Subset(slaveMeshID);
  if (slaveMeshAdapter != NULL) {
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;
    boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> slaveOperator = boost::dynamic_pointer_cast<
        AMP::Operator::MechanicsLinearFEOperator>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter,
                                                                                                 "MechanicsLinearFEOperator",
                                                                                                 input_db,
                                                                                                 slaveElementPhysicsModel));
    columnOperator->append(slaveOperator);

    boost::shared_ptr<AMP::Database> slaveSolver_db = columnPreconditioner_db->getDatabase("SlaveSolver"); 
    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> slaveSolverParams(new
        AMP::Solver::PetscKrylovSolverParameters(slaveSolver_db));
    slaveSolverParams->d_pOperator = slaveOperator;
//    slaveSolverParams->d_comm = globalComm;
    slaveSolverParams->d_comm = slaveMeshAdapter->getComm();
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> slaveSolver(new AMP::Solver::PetscKrylovSolver(slaveSolverParams));
    columnPreconditioner->append(slaveSolver);

    slaveLoadOperator = boost::dynamic_pointer_cast<
        AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter, 
                                                                                                 "SlaveLoadOperator", 
                                                                                                 input_db, 
                                                                                                 slaveElementPhysicsModel));
    AMP::LinearAlgebra::Variable::shared_ptr slaveVar = slaveOperator->getOutputVariable();
    slaveLoadOperator->setVariable(slaveVar);
  } // end if

  columnOperator->append(contactOperator);

  boost::shared_ptr<AMP::Database> contactPreconditioner_db = columnPreconditioner_db->getDatabase("ContactPreconditioner"); 
  boost::shared_ptr<AMP::Solver::MPCSolverParameters> contactPreconditionerParams(new 
      AMP::Solver::MPCSolverParameters(contactPreconditioner_db));
  contactPreconditionerParams->d_pOperator = contactOperator;
  boost::shared_ptr<AMP::Solver::MPCSolver> contactPreconditioner(new AMP::Solver::MPCSolver(contactPreconditionerParams));
  columnPreconditioner->append(contactPreconditioner);

  // Build a matrix shell operator to use the column operator with the petsc krylov solvers
  boost::shared_ptr<AMP::Database> matrixShellDatabase = input_db->getDatabase("MatrixShellOperator");
  boost::shared_ptr<AMP::Operator::OperatorParameters> matrixShellParams(new
      AMP::Operator::OperatorParameters(matrixShellDatabase));
  boost::shared_ptr<AMP::Operator::PetscMatrixShellOperator> matrixShellOperator(new
      AMP::Operator::PetscMatrixShellOperator(matrixShellParams));

  int numMasterLocalNodes = 0;
  int numSlaveLocalNodes = 0;
  if (masterMeshAdapter != NULL) { numMasterLocalNodes = masterMeshAdapter->numLocalElements(AMP::Mesh::Vertex); }
  if (slaveMeshAdapter != NULL) { numSlaveLocalNodes = slaveMeshAdapter->numLocalElements(AMP::Mesh::Vertex); }
  int matLocalSize = dofsPerNode * (numMasterLocalNodes + numSlaveLocalNodes);
  AMP_ASSERT( matLocalSize == dofManager->numLocalDOF() );
  if(!rank) {
    std::cout<<"numMasterNodes = "<<numMasterLocalNodes<<std::endl;
    std::cout<<"numSlaveNodes = "<<numSlaveLocalNodes<<std::endl;
    std::cout<<"matLocalSize = "<<matLocalSize<<std::endl;
  }
  matrixShellOperator->setComm(globalComm);
  matrixShellOperator->setMatLocalRowSize(matLocalSize);
  matrixShellOperator->setMatLocalColumnSize(matLocalSize);
  matrixShellOperator->setOperator(columnOperator); 

  AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = createVector(dofManager, columnVar, split);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = createVector(dofManager, columnVar, split);
  AMP::LinearAlgebra::Vector::shared_ptr columnResVec = createVector(dofManager, columnVar, split);
  columnSolVec->zero();
  columnRhsVec->zero();

  if (slaveLoadOperator != NULL) { slaveLoadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0); }


  columnSolVec->setToScalar(1.0);
  columnOperator->apply(nullVec, columnSolVec, columnResVec, 1.0, 0.0);
  contactOperator->addSlaveToMaster(columnResVec);
  contactOperator->setSlaveToZero(columnResVec);
  columnResVec->subtract(columnRhsVec, columnResVec);

  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = matrixShellOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = columnPreconditioner;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));
  linearSolver->setZeroInitialGuess(true);

  linearSolver->solve(columnRhsVec, columnSolVec);


#ifdef USE_SILO
  siloWriter->registerVector(columnSolVec, meshAdapter, AMP::Mesh::Vertex, "Solution");
  char outFileName[256];
  sprintf(outFileName, "MPC_%d", 0);
  siloWriter->writeFile(outFileName, 0);
#endif
  fout.close();

  ut->passes(exeName);
}

void myTest2(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName; 

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

//  int npes = globalComm.getSize();
  int rank = globalComm.getRank();

  // Load the input file
  globalComm.barrier();
  double inpReadBeginTime = MPI_Wtime();

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  globalComm.barrier();
  double inpReadEndTime = MPI_Wtime();
  if(!rank) {
    std::cout<<"Finished parsing the input file in "<<(inpReadEndTime - inpReadBeginTime)<<" seconds."<<std::endl;
  }

  // Load the meshes
  globalComm.barrier();
  double meshBeginTime = MPI_Wtime();

  AMP_INSIST(input_db->keyExists("FusedMesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("FusedMesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(globalComm);
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);

  globalComm.barrier();
  double meshEndTime = MPI_Wtime();
  if(!rank) {
    std::cout<<"Finished reading the mesh in "<<(meshEndTime - meshBeginTime)<<" seconds."<<std::endl;
  }

  
  int dofsPerNode = 3;
  int nodalGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr dofManager = AMP::Discretization::simpleDOFManager::create(meshAdapter,
      AMP::Mesh::Vertex, nodalGhostWidth, dofsPerNode, split);


  // build a column operator and a column preconditioner
  boost::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator(new AMP::Operator::ColumnOperator(emptyParams));

  boost::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase("LinearSolver"); 
  boost::shared_ptr<AMP::Database> columnPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::ColumnSolverParameters> columnPreconditionerParams(new
      AMP::Solver::ColumnSolverParameters(columnPreconditioner_db));
  columnPreconditionerParams->d_pOperator = columnOperator;
  boost::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

  // build the master and slave operators
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> masterOperator = boost::dynamic_pointer_cast<
        AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
                                                                                         "MasterBVPOperator",
                                                                                         input_db,
                                                                                         masterElementPhysicsModel));
    columnOperator->append(masterOperator);

    boost::shared_ptr<AMP::Database> masterSolver_db = columnPreconditioner_db->getDatabase("MasterSolver"); 
    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> masterSolverParams(new
        AMP::Solver::PetscKrylovSolverParameters(masterSolver_db));
    masterSolverParams->d_pOperator = masterOperator;
    masterSolverParams->d_comm = globalComm;
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> masterSolver(new AMP::Solver::PetscKrylovSolver(masterSolverParams));
    columnPreconditioner->append(masterSolver);

    boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> slaveLoadOperator = boost::dynamic_pointer_cast<
        AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, 
                                                                                                 "SlaveLoadOperator", 
                                                                                                 input_db, 
                                                                                                 masterElementPhysicsModel));
    AMP::LinearAlgebra::Variable::shared_ptr slaveVar = masterOperator->getOutputVariable();
    slaveLoadOperator->setVariable(slaveVar);


  // Build a matrix shell operator to use the column operator with the petsc krylov solvers
  boost::shared_ptr<AMP::Database> matrixShellDatabase = input_db->getDatabase("MatrixShellOperator");
  boost::shared_ptr<AMP::Operator::OperatorParameters> matrixShellParams(new
      AMP::Operator::OperatorParameters(matrixShellDatabase));
  boost::shared_ptr<AMP::Operator::PetscMatrixShellOperator> matrixShellOperator(new
      AMP::Operator::PetscMatrixShellOperator(matrixShellParams));

  int matLocalSize = dofsPerNode * meshAdapter->numLocalElements(AMP::Mesh::Vertex); 
  AMP_ASSERT( matLocalSize == dofManager->numLocalDOF() );
  matrixShellOperator->setComm(globalComm);
  matrixShellOperator->setMatLocalRowSize(matLocalSize);
  matrixShellOperator->setMatLocalColumnSize(matLocalSize);
  matrixShellOperator->setOperator(columnOperator); 

  AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = createVector(dofManager, columnVar, split);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = createVector(dofManager, columnVar, split);
  columnSolVec->zero();
  columnRhsVec->zero();

  slaveLoadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);

  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = matrixShellOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = columnPreconditioner;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));
  linearSolver->setZeroInitialGuess(true);

  linearSolver->solve(columnRhsVec, columnSolVec);


  ut->passes(exeName);
}






int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  AMP::UnitTest ut;

  std::vector<std::string> exeNames; 
  exeNames.push_back("testNodeToSegmentConstraintsOperator-cube");
//  exeNames.push_back("testNodeToSegmentConstraintsOperator-cylinder");
//  exeNames.push_back("testNodeToSegmentConstraintsOperator-pellet");

  try {
    for (size_t i = 0; i < exeNames.size(); ++i) { myTest(&ut, exeNames[i]); myTest2(&ut, exeNames[i]); }
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



