
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

  AMP_INSIST(input_db->keyExists("FusedMesh"), "Key ''FusedMesh'' is missing!");
  boost::shared_ptr<AMP::Database> fusedMesh_db = input_db->getDatabase("FusedMesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> fusedMeshParams(new AMP::Mesh::MeshParameters(fusedMesh_db));
  fusedMeshParams->setComm(globalComm);
  AMP::Mesh::Mesh::shared_ptr fusedMeshAdapter = AMP::Mesh::Mesh::buildMesh(fusedMeshParams);

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

  AMP::Discretization::DOFManager::shared_ptr fusedMeshDOFManager = AMP::Discretization::simpleDOFManager::create(fusedMeshAdapter,
      AMP::Mesh::Vertex, nodalGhostWidth, dofsPerNode, split);

/*  nodeToSegmentConstraintsOperator->reset(nodeToSegmentConstraintsOperatorParams);

  AMP::LinearAlgebra::Variable::shared_ptr dummyVariable(new AMP::LinearAlgebra::Variable("Dummy"));
  AMP::LinearAlgebra::Vector::shared_ptr dummyInVector = createVector(dofManager, dummyVariable, split);
  AMP::LinearAlgebra::Vector::shared_ptr dummyOutVector = createVector(dofManager, dummyVariable, split);

  nodeToSegmentConstraintsOperator->apply(dummyInVector, dummyInVector, dummyOutVector);
  nodeToSegmentConstraintsOperator->applyTranspose(dummyInVector, dummyInVector, dummyOutVector);*/


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

//  columnOperator->append(contactOperator);

  boost::shared_ptr<AMP::Database> contactPreconditioner_db = columnPreconditioner_db->getDatabase("ContactPreconditioner"); 
  boost::shared_ptr<AMP::Solver::MPCSolverParameters> contactPreconditionerParams(new 
      AMP::Solver::MPCSolverParameters(contactPreconditioner_db));
  contactPreconditionerParams->d_pOperator = contactOperator;
  boost::shared_ptr<AMP::Solver::MPCSolver> contactPreconditioner(new AMP::Solver::MPCSolver(contactPreconditionerParams));
  columnPreconditioner->append(contactPreconditioner);

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> fusedMeshElementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> fusedMeshOperator = boost::dynamic_pointer_cast<
      AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(fusedMeshAdapter,
                                                                                       "MasterBVPOperator",
                                                                                       input_db,
                                                                                       fusedMeshElementPhysicsModel));
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection>  fusedMeshLoadOperator = boost::dynamic_pointer_cast<
      AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(fusedMeshAdapter, 
                                                                                               "SlaveLoadOperator", 
                                                                                               input_db, 
                                                                                               fusedMeshElementPhysicsModel));
  AMP::LinearAlgebra::Variable::shared_ptr fusedMeshVar = fusedMeshOperator->getOutputVariable();
  fusedMeshLoadOperator->setVariable(fusedMeshVar);

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


  int fusedMeshMatLocalSize = dofsPerNode * fusedMeshAdapter->numLocalElements(AMP::Mesh::Vertex);
  AMP_ASSERT( fusedMeshMatLocalSize == fusedMeshDOFManager->numLocalDOF() );
  boost::shared_ptr<AMP::Operator::PetscMatrixShellOperator> fusedMeshMatrixShellOperator(new
      AMP::Operator::PetscMatrixShellOperator(matrixShellParams));
  fusedMeshMatrixShellOperator->setComm(globalComm);
  fusedMeshMatrixShellOperator->setMatLocalRowSize(fusedMeshMatLocalSize);
  fusedMeshMatrixShellOperator->setMatLocalColumnSize(fusedMeshMatLocalSize);
  fusedMeshMatrixShellOperator->setOperator(fusedMeshOperator);

  std::cout<<"MPC > NUMBER OF VERTICES IS "<<meshAdapter->numGlobalElements(AMP::Mesh::Vertex)<<std::endl;
  std::cout<<"MPC > NUMBER OF ELEMENTS IS "<<meshAdapter->numGlobalElements(AMP::Mesh::Volume)<<std::endl;
  std::cout<<"MPC > NUMBER OF CONSTRAINTS IS "<<contactOperator->numGlobalConstraints()<<std::endl;
  std::cout<<"FUSED > NUMBER OF VERTICES IS "<<fusedMeshAdapter->numGlobalElements(AMP::Mesh::Vertex)<<std::endl;
  std::cout<<"FUSED > NUMBER OF ELEMENTS IS "<<fusedMeshAdapter->numGlobalElements(AMP::Mesh::Volume)<<std::endl;


  AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = createVector(dofManager, columnVar, split);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = createVector(dofManager, columnVar, split);
  AMP::LinearAlgebra::Vector::shared_ptr columnResVec = createVector(dofManager, columnVar, split);
  columnSolVec->zero();
  columnRhsVec->zero();

  AMP::LinearAlgebra::Vector::shared_ptr fusedMeshSolVec = createVector(fusedMeshDOFManager, fusedMeshVar, split);
  AMP::LinearAlgebra::Vector::shared_ptr fusedMeshRhsVec = createVector(fusedMeshDOFManager, fusedMeshVar, split);
  AMP::LinearAlgebra::Vector::shared_ptr fusedMeshResVec = createVector(fusedMeshDOFManager, fusedMeshVar, split);
  fusedMeshSolVec->zero();
  fusedMeshRhsVec->zero();

  if (slaveLoadOperator != NULL) { slaveLoadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0); }

  fusedMeshLoadOperator->apply(nullVec, nullVec, fusedMeshRhsVec, 1.0, 0.0);

  std::cout<<std::setprecision(15);
  std::cout<<"MPC > RHS NORM = "<<columnRhsVec->L2Norm()<<std::endl;
  std::cout<<"FUSED > RHS NORM = "<<fusedMeshRhsVec->L2Norm()<<std::endl;

/*  columnSolVec->setToScalar(1.0);
  contactOperator->debugSetSlaveToZero(columnSolVec);
  std::cout<<"MPC > ZEROING SLAVE USING CONTACT OP PRIOR APPLY SOL L2NORM IS "<<columnSolVec->L2Norm()<<std::endl;
  columnSolVec->setToScalar(1.0);
  if (slaveMeshAdapter != NULL) {
    AMP::Mesh::MeshIterator slaveBoundaryIDIterator = slaveMeshAdapter->getBoundaryIDIterator(AMP::Mesh::Vertex, contactOperatorParams->d_SlaveBoundaryID);
    AMP::Mesh::MeshIterator boundaryIterator = slaveBoundaryIDIterator.begin(),
        boundaryIterator_end = slaveBoundaryIDIterator.end();
    for ( ; boundaryIterator != boundaryIterator_end; ++boundaryIterator) {
      std::vector<size_t> dofs;
      dofManager->getDOFs(boundaryIterator->globalID(), dofs);
      AMP_ASSERT( dofs.size() == dofsPerNode );
      std::vector<double> zeros(dofsPerNode, 0.0);
      columnSolVec->setLocalValuesByGlobalID(dofsPerNode, &(dofs[0]), &(zeros[0]));
    } // end for
  } // end if
  std::cout<<"MPC > ZEROING SLAVE USING MESH ITERATOR PRIOR APPLY SOL L2NORM IS "<<columnSolVec->L2Norm()<<std::endl;
  fusedMeshSolVec->setToScalar(1.0);
  std::cout<<"FUSED > PRIOR APPLY SOL L2NORM IS "<<fusedMeshSolVec->L2Norm()<<std::endl;*/

  AMP::LinearAlgebra::Vector::shared_ptr oldSolVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr oldResVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr dirVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr oldDirVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr matVec = columnSolVec->cloneVector();

  AMP::LinearAlgebra::Vector::shared_ptr fusedMeshOldSolVec = fusedMeshSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr fusedMeshOldResVec = fusedMeshSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr fusedMeshDirVec = fusedMeshSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr fusedMeshOldDirVec = fusedMeshSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr fusedMeshMatVec = fusedMeshSolVec->cloneVector();

  columnSolVec->setToScalar(1.0);
  columnOperator->apply(nullVec, columnSolVec, columnResVec, 1.0, 0.0);
  contactOperator->addSlaveToMaster(columnResVec);
  contactOperator->setSlaveToZero(columnResVec);
  columnResVec->subtract(columnRhsVec, columnResVec);

  fusedMeshSolVec->setToScalar(1.0);
  fusedMeshOperator->apply(nullVec, fusedMeshSolVec, fusedMeshResVec, 1.0, 0.0);
  fusedMeshResVec->subtract(fusedMeshRhsVec, fusedMeshResVec);

  oldDirVec->copyVector(columnResVec);
  fusedMeshOldDirVec->copyVector(fusedMeshResVec);

  int maxIters = input_db->getInteger("maxIters");
  for (int iter = 0; iter < maxIters; ++iter) {
    double resNorm = columnResVec->L2Norm();
    contactOperator->setSlaveToZero(columnSolVec);
    double solNorm = columnSolVec->L2Norm();
    contactOperator->copyMasterToSlave(columnSolVec);
    double fusedMeshResNorm = fusedMeshResVec->L2Norm();
    double fusedMeshSolNorm = fusedMeshSolVec->L2Norm();
    std::cout<<"Iter = "<<iter<<"  "
        <<"resNorm = "<<resNorm<<"  solNorm = "<<solNorm<<"\n"
        <<"fusedMeshResNorm = "<<fusedMeshResNorm<<"  fusedMeshSolNorm = "<<fusedMeshSolNorm<<"\n";

    contactOperator->copyMasterToSlave(oldDirVec);

    oldResVec->copyVector(columnResVec);
    oldSolVec->copyVector(columnSolVec);
    fusedMeshOldResVec->copyVector(fusedMeshResVec);
    fusedMeshOldSolVec->copyVector(fusedMeshSolVec);

    double alphaNumerator = oldResVec->dot(oldResVec);
    double fusedMeshAlphaNumerator = fusedMeshOldResVec->dot(fusedMeshOldResVec);

    columnOperator->apply(nullVec, oldDirVec, matVec, 1.0, 0.0);
    contactOperator->addSlaveToMaster(matVec);
    contactOperator->setSlaveToZero(matVec);
    fusedMeshOperator->apply(nullVec, fusedMeshOldDirVec, fusedMeshMatVec, 1.0, 0.0);

    double matVecNorm = matVec->L2Norm();
    double fusedMeshMatVecNorm = fusedMeshMatVec->L2Norm();
    std::cout<<"  matVecNorm = "<<matVecNorm<<"  "
      "fusedMeshMatVecNorm = "<<fusedMeshMatVecNorm<<"\n";

    double alphaDenominator = matVec->dot(oldDirVec);
    double fusedMeshAlphaDenominator = fusedMeshMatVec->dot(fusedMeshOldDirVec);

    double alpha = alphaNumerator / alphaDenominator;
    double fusedMeshAlpha = fusedMeshAlphaNumerator / fusedMeshAlphaDenominator;

    columnSolVec->axpy(alpha, oldDirVec, oldSolVec);
    columnResVec->axpy(-alpha, matVec, oldResVec);
    fusedMeshSolVec->axpy(fusedMeshAlpha, fusedMeshOldDirVec, fusedMeshOldSolVec);
    fusedMeshResVec->axpy(-fusedMeshAlpha, fusedMeshMatVec, fusedMeshOldResVec);

    double betaNumerator = columnResVec->dot(columnResVec);
    double fusedMeshBetaNumerator = fusedMeshResVec->dot(fusedMeshResVec);
    double beta = betaNumerator / alphaNumerator;
    double fusedMeshBeta = fusedMeshBetaNumerator / fusedMeshAlphaNumerator;
    std::cout<<"  beta = "<<beta<<"  "
        <<"fusedMeshBeta = "<<fusedMeshBeta<<"\n";

    dirVec->axpy(beta, oldDirVec, columnResVec);
    fusedMeshDirVec->axpy(fusedMeshBeta, fusedMeshOldDirVec, fusedMeshResVec);

    oldDirVec->copyVector(dirVec);
    fusedMeshOldDirVec->copyVector(fusedMeshDirVec);
  } // end for

/*
  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = matrixShellOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = columnPreconditioner;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));
  linearSolver->setZeroInitialGuess(true);

  std::cout<<"MPC > "<<std::endl;
  linearSolver->solve(columnRhsVec, columnSolVec);


  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> fusedMeshLinearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));
  fusedMeshLinearSolver->setZeroInitialGuess(true);

  std::cout<<"FUSED > "<<std::endl;
  fusedMeshLinearSolver->solve(fusedMeshRhsVec, fusedMeshSolVec);
*/
#ifdef USE_SILO
  siloWriter->registerVector(columnSolVec, meshAdapter, AMP::Mesh::Vertex, "Solution");
  char outFileName[256];
  sprintf(outFileName, "MPC_%d", 0);
  siloWriter->writeFile(outFileName, 0);
#endif
  fout.close();

  ut->passes(exeName);
}






int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
//  boost::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit( new AMP::Mesh::initializeLibMesh(globalComm) );
  AMP::UnitTest ut;

  std::vector<std::string> exeNames; 
  exeNames.push_back("testNodeToSegmentConstraintsOperator-cube");
  exeNames.push_back("testNodeToSegmentConstraintsOperator-cylinder");
  exeNames.push_back("testNodeToSegmentConstraintsOperator-pellet");

  try {
    for (size_t i = 0; i < exeNames.size(); ++i) { myTest(&ut, exeNames[i]); }
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  ut.report();
  int num_failed = ut.NumFailGlobal();

//  libmeshInit.reset();
  AMP::AMPManager::shutdown();
  return num_failed;
}  



