
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
#include "ampmesh/latex_visualization_tools.h"
#include "ampmesh/euclidean_geometry_tools.h"

void drawBoundaryID(AMP::Mesh::Mesh::shared_ptr meshAdapter, int boundaryID, std::ostream &os, double const * point_of_view) {
  AMP::Mesh::MeshIterator boundaryIterator = meshAdapter->getBoundaryIDIterator(AMP::Mesh::Face, boundaryID);
  AMP::Mesh::MeshIterator boundaryIterator_begin = boundaryIterator.begin(), 
      boundaryIterator_end = boundaryIterator.end();
  std::vector<AMP::Mesh::MeshElement> faceVertices;
  std::vector<double> faceVertexCoordinates;
  double faceData[12];
  double const * faceDataPtr[4] = { faceData, faceData+3, faceData+6, faceData+9 };

  os<<std::setprecision(6)<<std::fixed;

  for (boundaryIterator = boundaryIterator_begin; boundaryIterator != boundaryIterator_end; ++boundaryIterator) {
    faceVertices = boundaryIterator->getElements(AMP::Mesh::Vertex);
    AMP_ASSERT( faceVertices.size() == 4 );
    for (size_t i = 0; i < 4; ++i) {
      faceVertexCoordinates = faceVertices[i].coord();
      AMP_ASSERT( faceVertexCoordinates.size() == 3 );
      std::copy(faceVertexCoordinates.begin(), faceVertexCoordinates.end(), faceData+3*i);
    } // end for i
    triangle_t t(faceDataPtr[0], faceDataPtr[1], faceDataPtr[2]);

    std::string option = "";
    if (compute_scalar_product(point_of_view, t.get_normal()) > 0.0) {
      os<<"\\draw["<<option<<"]\n";
      write_face(faceDataPtr, os);
    } // end if
  } // end for
}

void myPCG(AMP::LinearAlgebra::Vector::shared_ptr rhs, AMP::LinearAlgebra::Vector::shared_ptr sol, 
    AMP::Operator::Operator::shared_ptr op, boost::shared_ptr<AMP::Solver::SolverStrategy> pre,
    size_t maxIters, double tol, bool verbose = false, std::ostream &os = std::cout) {
  AMP::LinearAlgebra::Vector::shared_ptr res = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr dir = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr ext = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr oldSol = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr oldRes = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr oldDir = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr matVec = sol->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  op->apply(nullVec, sol, matVec, 1.0, 0.0);
  oldRes->subtract(rhs, matVec);
  pre->solve(oldRes, ext);
  oldDir->copyVector(ext);
  oldSol->copyVector(sol);
  double initialResNorm = oldRes->L2Norm();
  if (verbose) { os<<std::setprecision(15)<<"  iter=0  itialResNorm="<<initialResNorm<<"\n"; }
  for (size_t iter = 0; iter < maxIters; ++iter) {
    if (verbose) { os<<"  iter="<<iter+1<<"  "; }
    op->apply(nullVec, oldDir, matVec, 1.0, 0.0);
    double extDOToldRes = ext->dot(oldRes);
    double oldDirDOTmatVec = oldDir->dot(matVec);
    double alpha = extDOToldRes / oldDirDOTmatVec;
    if (verbose) { os<<"alpha="<<alpha<<"  "; }
    if (verbose) { os<<"oldDirDOTmatVec="<<oldDirDOTmatVec<<"  "; }
    sol->axpy(alpha, oldDir, oldSol);
    res->axpy(-alpha, matVec, oldRes);
    double resNorm = res->L2Norm();
    if (verbose) { os<<"resNorm="<<resNorm<<"  "; }
    pre->solve(res, ext);
    double extDOTres = ext->dot(res);
    double beta = extDOTres / extDOToldRes;
    if (verbose) { os<<"beta="<<beta<<"\n"; }
    dir->axpy(beta, oldDir, ext);
    oldSol->copyVector(sol);
    oldRes->copyVector(res);
    oldDir->copyVector(dir);
  } // end for
}


void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName; 

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

#ifdef USE_SILO
  AMP::Mesh::SiloIO::shared_ptr siloWriter(new AMP::Mesh::SiloIO);
#endif

//  int npes = globalComm.getSize();
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

  // Create a DOF manager
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

  // Build the contact operator
  AMP_INSIST(input_db->keyExists("ContactOperator"), "Key ''ContactOperator'' is missing!");
  boost::shared_ptr<AMP::Database> contact_db = input_db->getDatabase("ContactOperator");
  boost::shared_ptr<AMP::Operator::NodeToSegmentConstraintsOperatorParameters> 
      contactOperatorParams( new AMP::Operator::NodeToSegmentConstraintsOperatorParameters(contact_db) );
  contactOperatorParams->d_DOFsPerNode = dofsPerNode;
  contactOperatorParams->d_DOFManager = dofManager;
  contactOperatorParams->d_GlobalComm = globalComm;
  contactOperatorParams->d_Mesh = meshAdapter;

  boost::shared_ptr<AMP::Operator::NodeToSegmentConstraintsOperator> 
      contactOperator( new AMP::Operator::NodeToSegmentConstraintsOperator(contactOperatorParams) );

  // TODO: RESET IN CONSTRUCTOR?
  contactOperator->reset(contactOperatorParams);

  // build the master and slave operators
  AMP::Mesh::MeshID masterMeshID = contactOperator->getMasterMeshID();
  AMP::Mesh::Mesh::shared_ptr masterMeshAdapter = meshAdapter->Subset(masterMeshID);
  if (masterMeshAdapter.get() != NULL) {
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
double point_of_view[3] = { 1.0, 1.0, 1.0 };
drawBoundaryID(masterMeshAdapter, 1, fout, point_of_view);
drawBoundaryID(masterMeshAdapter, 4, fout, point_of_view);
  } // end if

  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> slaveLoadOperator;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> slaveBVPOperator;

  AMP::Mesh::MeshID slaveMeshID = contactOperator->getSlaveMeshID();
  AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter = meshAdapter->Subset(slaveMeshID);
  if (slaveMeshAdapter.get() != NULL) {
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;

    boost::shared_ptr<AMP::Database> slaveSolver_db = columnPreconditioner_db->getDatabase("SlaveSolver"); 
    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> slaveSolverParams(new
        AMP::Solver::PetscKrylovSolverParameters(slaveSolver_db));

    bool useSlaveBVPOperator = input_db->getBool("useSlaveBVPOperator");
    if (useSlaveBVPOperator) {
      slaveBVPOperator = boost::dynamic_pointer_cast<
          AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter,
                                                                                           "SlaveBVPOperator",
                                                                                           input_db,
                                                                                           slaveElementPhysicsModel));
      columnOperator->append(slaveBVPOperator);
      slaveSolverParams->d_pOperator = slaveBVPOperator;
    } else {
      boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> slaveMechanicsLinearFEOperator = boost::dynamic_pointer_cast<
          AMP::Operator::MechanicsLinearFEOperator>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter,
                                                                                                   "MechanicsLinearFEOperator",
                                                                                                   input_db,
                                                                                                   slaveElementPhysicsModel));
      columnOperator->append(slaveMechanicsLinearFEOperator);
 
      slaveSolverParams->d_pOperator = slaveMechanicsLinearFEOperator;

      slaveLoadOperator = boost::dynamic_pointer_cast<
          AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter, 
                                                                                                   "SlaveLoadOperator", 
                                                                                                   input_db, 
                                                                                                   slaveElementPhysicsModel));
      AMP::LinearAlgebra::Variable::shared_ptr slaveVar = slaveMechanicsLinearFEOperator->getOutputVariable();
      slaveLoadOperator->setVariable(slaveVar);
    } // end if

//    slaveSolverParams->d_comm = globalComm;
    slaveSolverParams->d_comm = slaveMeshAdapter->getComm();
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> slaveSolver(new AMP::Solver::PetscKrylovSolver(slaveSolverParams));
    columnPreconditioner->append(slaveSolver);

  } // end if

  boost::shared_ptr<AMP::Database> contactPreconditioner_db = columnPreconditioner_db->getDatabase("ContactPreconditioner"); 
  boost::shared_ptr<AMP::Solver::MPCSolverParameters> contactPreconditionerParams(new 
      AMP::Solver::MPCSolverParameters(contactPreconditioner_db));
  contactPreconditionerParams->d_pOperator = contactOperator;
  boost::shared_ptr<AMP::Solver::MPCSolver> contactPreconditioner(new AMP::Solver::MPCSolver(contactPreconditionerParams));
  columnPreconditioner->append(contactPreconditioner);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = createVector(dofManager, columnVar, split);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = createVector(dofManager, columnVar, split);
  columnSolVec->zero();
  columnRhsVec->zero();

  // compute f
  if (slaveLoadOperator.get() != NULL) { 
    slaveLoadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);
  } // end if

  // get d
  AMP::LinearAlgebra::Vector::shared_ptr tmpVec = createVector(dofManager, columnVar, split);
  tmpVec->zero();
  contactOperator->addShiftToSlave(tmpVec);
  double tmpVecRes = tmpVec->L2Norm();
  if(!rank) std::cout<<"tmpVecRes="<<tmpVecRes<<"\n";

  // compute - Kd
  AMP::LinearAlgebra::Vector::shared_ptr rhsCorrectionVec = createVector(dofManager, columnVar, split);
  columnOperator->apply(nullVec, tmpVec, rhsCorrectionVec, -1.0, 0.0);
  columnOperator->append(contactOperator);
  double rhsCorrectionVecNorm = rhsCorrectionVec->L2Norm();
  if(!rank) std::cout<<"rhsCorrectionVecNorm="<<rhsCorrectionVecNorm<<"\n";

  // f = f - Kd
  columnRhsVec->add(columnRhsVec, rhsCorrectionVec);

  // f^m = f^m + C^T f^s
  // f^s = 0
  contactOperator->addSlaveToMaster(columnRhsVec);
  contactOperator->setSlaveToZero(columnRhsVec);

  // apply dirichlet rhs correction
  if (slaveBVPOperator.get() != NULL) {
    slaveBVPOperator->modifyRHSvector(columnRhsVec);
  } // end if


  bool usePetscKrylovSolver = input_db->getBool("usePetscKrylovSolver");
  if (usePetscKrylovSolver) {
  // Build a matrix shell operator to use the column operator with the petsc krylov solvers
  boost::shared_ptr<AMP::Database> matrixShellDatabase = input_db->getDatabase("MatrixShellOperator");
  boost::shared_ptr<AMP::Operator::OperatorParameters> matrixShellParams(new
      AMP::Operator::OperatorParameters(matrixShellDatabase));
  boost::shared_ptr<AMP::Operator::PetscMatrixShellOperator> matrixShellOperator(new
      AMP::Operator::PetscMatrixShellOperator(matrixShellParams));

  int numMasterLocalNodes = 0;
  int numSlaveLocalNodes = 0;
  if (masterMeshAdapter.get() != NULL) { numMasterLocalNodes = masterMeshAdapter->numLocalElements(AMP::Mesh::Vertex); }
  if (slaveMeshAdapter.get() != NULL) { numSlaveLocalNodes = slaveMeshAdapter->numLocalElements(AMP::Mesh::Vertex); }
  int matLocalSize = dofsPerNode * (numMasterLocalNodes + numSlaveLocalNodes);
  AMP_ASSERT( matLocalSize == dofManager->numLocalDOF() );
  matrixShellOperator->setComm(globalComm);
  matrixShellOperator->setMatLocalRowSize(matLocalSize);
  matrixShellOperator->setMatLocalColumnSize(matLocalSize);
  matrixShellOperator->setOperator(columnOperator); 

  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = matrixShellOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = columnPreconditioner;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));
  linearSolver->setZeroInitialGuess(true);

  linearSolver->solve(columnRhsVec, columnSolVec);

  } else {
    size_t myPCGmaxIters = input_db->getInteger("myPCGmaxIters");
    myPCG(columnRhsVec, columnSolVec, columnOperator, columnPreconditioner, myPCGmaxIters, -99.0, true); 
  }
  // u^s = C u^m + d
  contactOperator->copyMasterToSlave(columnSolVec);
  contactOperator->addShiftToSlave(columnSolVec);

  meshAdapter->displaceMesh(columnSolVec);

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
  if (!rank) { std::cout<<"### FUSED MESHES CASE ###"<<std::endl; }

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


  AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = createVector(dofManager, columnVar, split);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = createVector(dofManager, columnVar, split);
  columnSolVec->zero();
  columnRhsVec->zero();

  slaveLoadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);

  bool usePetscKrylovSolver = input_db->getBool("usePetscKrylovSolver");
  if (usePetscKrylovSolver) {
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

  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = matrixShellOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = columnPreconditioner;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));
  linearSolver->setZeroInitialGuess(true);

  linearSolver->solve(columnRhsVec, columnSolVec);

  } else {
    size_t myPCGmaxIters = input_db->getInteger("myPCGmaxIters");
    myPCG(columnRhsVec, columnSolVec, columnOperator, columnPreconditioner, myPCGmaxIters, -99.0, true); 
  }

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
    for (size_t i = 0; i < exeNames.size(); ++i) { 
      myTest(&ut, exeNames[i]); 
//      myTest2(&ut, exeNames[i]); 
    } // end for
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



