
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
#include "utils/Writer.h"

#include "operators/OperatorBuilder.h"
#include "operators/LinearBVPOperator.h"
#include "operators/ColumnOperator.h"
#include "operators/petsc/PetscMatrixShellOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/mechanics/MechanicsModelParameters.h"
#include "operators/mechanics/MechanicsMaterialModel.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/contact/NodeToFaceContactOperator.h"
#include "operators/CustomConstraintsEliminationOperator.h"
#include "operators/mechanics/IsotropicElasticModel.h"

#include "solvers/ColumnSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/ConstraintsEliminationSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"


#include "utils/ReadTestMesh.h"

#include <fstream>
#include <boost/lexical_cast.hpp>
#include "ampmesh/latex_visualization_tools.h"
#include "ampmesh/euclidean_geometry_tools.h"

#include "testNodeToFaceContactOperator.h"

void selectNodes(AMP::Mesh::Mesh::shared_ptr mesh, std::vector<AMP::Mesh::MeshElementID> & nodesGlobalIDs) {
  AMP::Mesh::MeshIterator meshIterator = mesh->getBoundaryIDIterator(AMP::Mesh::Vertex, 3);
  AMP::Mesh::MeshIterator meshIterator_begin = meshIterator.begin();
  AMP::Mesh::MeshIterator meshIterator_end = meshIterator.end();
  nodesGlobalIDs.clear();
  for (meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator) {
    std::vector<double> coord = meshIterator->coord();
    if (std::abs(coord[2] - 0.5) < 1.0e-14) {
      nodesGlobalIDs.push_back(meshIterator->globalID());
//      std::cout<<nodesGlobalIDs.size()<<"  ("<<coord[0]<<", "<<coord[1]<<", "<<coord[2]<<")"<<std::endl;
    } // end if
  } // end for
}

void printNodesValues(AMP::Mesh::Mesh::shared_ptr mesh, std::vector<AMP::Mesh::MeshElementID> const & nodesGlobalIDs, AMP::LinearAlgebra::Vector::shared_ptr vectorField1, AMP::LinearAlgebra::Vector::shared_ptr vectorField2, std::ostream & os = std::cout) {
  AMP::Discretization::DOFManager::shared_ptr dofManager = vectorField1->getDOFManager();
  for (size_t i = 0; i < nodesGlobalIDs.size(); ++i) {
    std::vector<double> coord = mesh->getElement(nodesGlobalIDs[i]).coord();
    std::vector<size_t> dofIndices;
    dofManager->getDOFs(nodesGlobalIDs[i], dofIndices);
    AMP_ASSERT(dofIndices.size() == 1);
    double value1 = vectorField1->getLocalValueByGlobalID(dofIndices[0]);
    double value2 = vectorField2->getLocalValueByGlobalID(dofIndices[0]);
    os<<std::setprecision(15)<<coord[0]<<"  "<<value1<<"  "<<value2<<"\n";
  } // end for i 
}

void getConcentratedLoadAtNodes(double loadParameter, double loadCutoff, AMP::Mesh::Mesh::shared_ptr meshAdapter, AMP::LinearAlgebra::Vector::shared_ptr loadVector, AMP::Discretization::DOFManager::shared_ptr dofManager) {
  static std::vector<double> loadValues;
  static std::vector<size_t> dofIndices;
  AMP_ASSERT(loadValues.size() == dofIndices.size());

  if (loadValues.empty()) {
    double totalLoad = 0.0;
    AMP::Mesh::MeshIterator boundaryIterator = meshAdapter->getBoundaryIDIterator(AMP::Mesh::Vertex, 4, 0);
    AMP::Mesh::MeshIterator boundaryIterator_begin = boundaryIterator.begin(),
      boundaryIterator_end = boundaryIterator.end();
    std::vector<double> vertexCoordinates;
    std::vector<size_t> vertexDofIndices;
    for (boundaryIterator = boundaryIterator_begin; boundaryIterator != boundaryIterator_end; ++boundaryIterator) {
      vertexCoordinates = boundaryIterator->coord();
      AMP_ASSERT( vertexCoordinates.size() == 3 );
      dofManager->getDOFs(boundaryIterator->globalID(), vertexDofIndices);
      AMP_ASSERT( vertexDofIndices.size() == 3 );

      if (vertexCoordinates[1] > loadCutoff) {
        loadValues.push_back(loadParameter);
        dofIndices.push_back(vertexDofIndices[1]);
        if (vertexCoordinates[1] - 1.49213 < 0.0) {
          loadValues.back() /= 2.0; 
        } // end if
        if ((std::abs(vertexCoordinates[2] - 0.0) < 1.0e-14)
          || (std::abs(vertexCoordinates[2] - 1.0) < 1.0e-14)) {
          loadValues.back() /= 2.0; 
        } // end if
        totalLoad += loadValues.back();
//        std::cout<<loadValues.size()<<"  "<<loadValues.back()<<"  "<<dofIndices.back()<<" ("<<vertexCoordinates[0]<<", "<<vertexCoordinates[1]<<" ,"<<vertexCoordinates[2]<<")\n";
      } // end if
    } // end for
    std::cout<<"TOTAL load="<<totalLoad<<"\n";
    AMP_ASSERT(loadValues.size() > 0);
  } // end if

  loadVector->zero();
  loadVector->setLocalValuesByGlobalID(loadValues.size(), &(dofIndices[0]), &(loadValues[0]));
  loadVector->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
}

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName; 

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

#ifdef USE_EXT_SILO
  AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
  siloWriter->setDecomposition(1);
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
  AMP::Discretization::DOFManager::shared_ptr dispDofManager = AMP::Discretization::simpleDOFManager::create(meshAdapter,
      AMP::Mesh::Vertex, nodalGhostWidth, dofsPerNode, split);

  // Build a column operator and a column preconditioner
  boost::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator(new AMP::Operator::ColumnOperator(emptyParams));

  boost::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase("LinearSolver"); 
  boost::shared_ptr<AMP::Database> columnPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::ColumnSolverParameters> columnPreconditionerParams(new
      AMP::Solver::ColumnSolverParameters(columnPreconditioner_db));
  columnPreconditionerParams->d_pOperator = columnOperator;
  boost::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

  // Get the mechanics material model for the contact operator
  boost::shared_ptr<AMP::Database> model_db = input_db->getDatabase("MasterMechanicsMaterialModel");
  boost::shared_ptr<AMP::Operator::MechanicsModelParameters> masterMechanicsMaterialModelParams(new AMP::Operator::MechanicsModelParameters(model_db));
  boost::shared_ptr<AMP::Operator::MechanicsMaterialModel> masterMechanicsMaterialModel(new AMP::Operator::IsotropicElasticModel(masterMechanicsMaterialModelParams));

  // ... needed for computing stresses
  boost::shared_ptr<AMP::Database> slaveMechanicsMaterialModel_db = input_db->getDatabase("SlaveMechanicsMaterialModel");
  boost::shared_ptr<AMP::Operator::MechanicsModelParameters> slaveMechanicsMaterialModelParams(new AMP::Operator::MechanicsModelParameters(slaveMechanicsMaterialModel_db));
  boost::shared_ptr<AMP::Operator::MechanicsMaterialModel> slaveMechanicsMaterialModel(new AMP::Operator::IsotropicElasticModel(slaveMechanicsMaterialModelParams));

  // Build the contact operator
  AMP_INSIST(input_db->keyExists("ContactOperator"), "Key ''ContactOperator'' is missing!");
  boost::shared_ptr<AMP::Database> contact_db = input_db->getDatabase("ContactOperator");
  boost::shared_ptr<AMP::Operator::ContactOperatorParameters> 
      contactOperatorParams( new AMP::Operator::ContactOperatorParameters(contact_db) );
  contactOperatorParams->d_DOFsPerNode = dofsPerNode;
  contactOperatorParams->d_DOFManager = dispDofManager;
  contactOperatorParams->d_GlobalComm = globalComm;
  contactOperatorParams->d_Mesh = meshAdapter;
  contactOperatorParams->d_MasterMechanicsMaterialModel = masterMechanicsMaterialModel;
  contactOperatorParams->reset(); // got segfault at constructor since d_Mesh was pointing to NULL

  boost::shared_ptr<AMP::Operator::NodeToFaceContactOperator> 
      contactOperator( new AMP::Operator::NodeToFaceContactOperator(contactOperatorParams) );

  contactOperator->initialize();
  contactOperator->setContactIsFrictionless(contact_db->getBoolWithDefault("ContactIsFrictionless", false));
  
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> masterBVPOperator;

  // Build the master and slave operators
  AMP::Mesh::MeshID masterMeshID = contactOperator->getMasterMeshID();
  AMP::Mesh::Mesh::shared_ptr masterMeshAdapter = meshAdapter->Subset(masterMeshID);
  if (masterMeshAdapter.get() != NULL) {
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
    masterBVPOperator = boost::dynamic_pointer_cast<
        AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(masterMeshAdapter,
                                                                                         "MasterBVPOperator",
                                                                                         input_db,
                                                                                         masterElementPhysicsModel));
    columnOperator->append(masterBVPOperator);

    boost::shared_ptr<AMP::Database> masterSolver_db = columnPreconditioner_db->getDatabase("MasterSolver"); 
    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> masterSolverParams(new
        AMP::Solver::PetscKrylovSolverParameters(masterSolver_db));
    masterSolverParams->d_pOperator = masterBVPOperator;
    masterSolverParams->d_comm = masterMeshAdapter->getComm();
//    masterSolverParams->d_comm = globalComm;
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> masterSolver(new AMP::Solver::PetscKrylovSolver(masterSolverParams));
    columnPreconditioner->append(masterSolver);
//    boost::shared_ptr<AMP::Database> masterSolver_db = columnPreconditioner_db->getDatabase("NewMasterSolver"); 
//    boost::shared_ptr<AMP::Solver::SolverStrategyParameters> masterSolverParams(new AMP::Solver::SolverStrategyParameters(masterSolver_db));
//    masterSolverParams->d_pOperator = masterBVPOperator;
//    boost::shared_ptr<AMP::Solver::TrilinosMLSolver> masterSolver(new AMP::Solver::TrilinosMLSolver(masterSolverParams));
//    columnPreconditioner->append(masterSolver);


//std::fstream masterFout;
//masterFout.open("master_pellet", std::fstream::out);
//double point_of_view[3] = { 1.0, 1.0, 1.0 };
//drawFacesOnBoundaryID(masterMeshAdapter, 0, masterFout, point_of_view, "blue");
//drawFacesOnBoundaryID(masterMeshAdapter, 1, masterFout, point_of_view, "green");
//drawFacesOnBoundaryID(masterMeshAdapter, 2, masterFout, point_of_view, "red");
//drawFacesOnBoundaryID(masterMeshAdapter, 3, masterFout, point_of_view, "magenta");
//drawFacesOnBoundaryID(masterMeshAdapter, 4, masterFout, point_of_view, "black");
//drawFacesOnBoundaryID(masterMeshAdapter, 5, masterFout, point_of_view, "orange");
//drawFacesOnBoundaryID(masterMeshAdapter, 6, masterFout, point_of_view, "pink");
//drawFacesOnBoundaryID(masterMeshAdapter, 7, masterFout, point_of_view, "violet");
////drawFacesOnBoundaryID(masterMeshAdapter, 1, masterFout, point_of_view);
////drawFacesOnBoundaryID(masterMeshAdapter, 4, masterFout, point_of_view);
//masterFout.close();
  } // end if

  boost::shared_ptr<AMP::Operator::LinearBVPOperator> slaveBVPOperator;

  AMP::Mesh::MeshID slaveMeshID = contactOperator->getSlaveMeshID();
  AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter = meshAdapter->Subset(slaveMeshID);
  if (slaveMeshAdapter.get() != NULL) {
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;

    slaveBVPOperator = boost::dynamic_pointer_cast<
        AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter,
                                                                                         "SlaveBVPOperator",
                                                                                         input_db,
                                                                                         slaveElementPhysicsModel));
    columnOperator->append(slaveBVPOperator);

    boost::shared_ptr<AMP::Database> slaveSolver_db = columnPreconditioner_db->getDatabase("SlaveSolver"); 
    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> slaveSolverParams(new
        AMP::Solver::PetscKrylovSolverParameters(slaveSolver_db));
    slaveSolverParams->d_pOperator = slaveBVPOperator;
    slaveSolverParams->d_comm = slaveMeshAdapter->getComm();
//    slaveSolverParams->d_comm = globalComm;
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> slaveSolver(new AMP::Solver::PetscKrylovSolver(slaveSolverParams));
    columnPreconditioner->append(slaveSolver);
//    boost::shared_ptr<AMP::Database> slaveSolver_db = columnPreconditioner_db->getDatabase("NewSlaveSolver"); 
//    boost::shared_ptr<AMP::Solver::SolverStrategyParameters> slaveSolverParams(new AMP::Solver::SolverStrategyParameters(slaveSolver_db));
//    slaveSolverParams->d_pOperator = slaveBVPOperator;
//    boost::shared_ptr<AMP::Solver::TrilinosMLSolver> slaveSolver(new AMP::Solver::TrilinosMLSolver(slaveSolverParams));
//    columnPreconditioner->append(slaveSolver);

//std::fstream slaveFout;
//slaveFout.open("slave_pellet", std::fstream::out);
//double point_of_view[3] = { 1.0, 1.0, 1.0 };
//drawFacesOnBoundaryID(slaveMeshAdapter, 0, slaveFout, point_of_view, "dashed,red");
////drawFacesOnBoundaryID(slaveMeshAdapter, 1, slaveFout, point_of_view, "dashed");
////drawFacesOnBoundaryID(slaveMeshAdapter, 4, slaveFout, point_of_view, "dashed");
////drawVerticesOnBoundaryID(slaveMeshAdapter, 2, slaveFout, point_of_view, "red");
//slaveFout.close();
  } // end if

  boost::shared_ptr<AMP::InputDatabase> dummyOperator_db(new AMP::InputDatabase("dummyOperator_db"));
  boost::shared_ptr<AMP::Operator::OperatorParmeters> dummyOperatorParams(new AMP::Operator::OperatorParmeters(dummyOperator_db));
  boost::shared_ptr<AMP::Operator::CustomConstraintsEliminationOperator> dirichletCustomConstraintsEliminationOperator(new AMP::Operator::CustomConstraintsEliminationOperator(dummyOperatorParams));

  boost::shared_ptr<AMP::Database> contactPreconditioner_db = columnPreconditioner_db->getDatabase("ContactPreconditioner"); 
  boost::shared_ptr<AMP::Solver::ConstraintsEliminationSolverParameters> contactPreconditionerParams(new 
      AMP::Solver::ConstraintsEliminationSolverParameters(contactPreconditioner_db));
  contactPreconditionerParams->d_pOperator = contactOperator;
  boost::shared_ptr<AMP::Solver::ConstraintsEliminationSolver> contactPreconditioner(new AMP::Solver::ConstraintsEliminationSolver(contactPreconditionerParams));
  columnPreconditioner->append(contactPreconditioner);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = createVector(dispDofManager, columnVar, split);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = createVector(dispDofManager, columnVar, split);
  columnSolVec->zero();
  columnRhsVec->zero();

  AMP::LinearAlgebra::Variable::shared_ptr tempVar(new AMP::LinearAlgebra::Variable("temperature"));
  AMP::LinearAlgebra::Variable::shared_ptr dispVar = columnOperator->getOutputVariable();
  AMP::Discretization::DOFManager::shared_ptr tempDofManager = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth, 1 , split);
  AMP::LinearAlgebra::Vector::shared_ptr tempVec = AMP::LinearAlgebra::createVector(tempDofManager, tempVar, split);
  double const referenceTemperature = 300.0;
  tempVec->setToScalar(referenceTemperature);
  double const thermalExpansionCoefficient = 2.0e-4;

  AMP::LinearAlgebra::Vector::shared_ptr sigma_xx = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_xx")), split);
  AMP::LinearAlgebra::Vector::shared_ptr sigma_yy = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_yy")), split);
  AMP::LinearAlgebra::Vector::shared_ptr sigma_zz = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_zz")), split);
  AMP::LinearAlgebra::Vector::shared_ptr sigma_yz = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_yz")), split);
  AMP::LinearAlgebra::Vector::shared_ptr sigma_xz = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_xz")), split);
  AMP::LinearAlgebra::Vector::shared_ptr sigma_xy = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_xy")), split);
  AMP::LinearAlgebra::Vector::shared_ptr sigma_eff = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_eff")), split);
  AMP::LinearAlgebra::Vector::shared_ptr activeSetVec = sigma_eff->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr suckItVec = sigma_eff->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr lickItVec = sigma_eff->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr surfaceTractionVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr normalVectorVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr contactShiftVec = columnSolVec->cloneVector();
  activeSetVec->zero();
  suckItVec->zero();
  lickItVec->zero();
  surfaceTractionVec->zero();
  normalVectorVec->zero();
  contactShiftVec->zero();

  computeStressTensor(meshAdapter, columnSolVec, 
      sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy,
      sigma_eff, 1.0e6, 0.3,
      referenceTemperature, thermalExpansionCoefficient, tempVec);

  bool skipDisplaceMesh = true;
  contactOperator->updateActiveSet(nullVec, skipDisplaceMesh);

  {//TODO: this is tmp only
  normalVectorVec->zero();
  std::vector<AMP::Mesh::MeshElementID> const * pointerToActiveSet;
  contactOperator->getActiveSet(pointerToActiveSet);
  std::vector<size_t> activeSetDispDOFsIndicesBeforeUpdate;
  dispDofManager->getDOFs(*pointerToActiveSet, activeSetDispDOFsIndicesBeforeUpdate);
  std::vector<double> const * slaveVerticesNormalVector;
  std::vector<double> const * slaveVerticesSurfaceTraction;
  contactOperator->getSlaveVerticesNormalVectorAndSurfaceTraction(slaveVerticesNormalVector, slaveVerticesSurfaceTraction);
//  AMP_ASSERT( slaveVerticesNormalVector->size() == 3*sizeOfActiveSetBeforeUpdate );
//  AMP_ASSERT( slaveVerticesSurfaceTraction->size() == 3*sizeOfActiveSetBeforeUpdate );
  normalVectorVec->setLocalValuesByGlobalID(slaveVerticesNormalVector->size(), &(activeSetDispDOFsIndicesBeforeUpdate[0]), &((*slaveVerticesNormalVector)[0]));
  }
#ifdef USE_EXT_SILO
  {
    siloWriter->registerVector(columnSolVec, meshAdapter, AMP::Mesh::Vertex, "Solution");
    siloWriter->registerVector(sigma_eff, meshAdapter, AMP::Mesh::Vertex, "vonMises");
    siloWriter->registerVector(sigma_xx, meshAdapter, AMP::Mesh::Vertex, "sigma_xx");
    siloWriter->registerVector(sigma_yy, meshAdapter, AMP::Mesh::Vertex, "sigma_yy");
    siloWriter->registerVector(sigma_zz, meshAdapter, AMP::Mesh::Vertex, "sigma_zz");
    siloWriter->registerVector(sigma_yz, meshAdapter, AMP::Mesh::Vertex, "sigma_yz");
    siloWriter->registerVector(sigma_xz, meshAdapter, AMP::Mesh::Vertex, "sigma_xz");
    siloWriter->registerVector(sigma_xy, meshAdapter, AMP::Mesh::Vertex, "sigma_xy");
    siloWriter->registerVector(activeSetVec, meshAdapter, AMP::Mesh::Vertex, "Contact");
    siloWriter->registerVector(surfaceTractionVec, meshAdapter, AMP::Mesh::Vertex, "Traction");
    siloWriter->registerVector(normalVectorVec, meshAdapter, AMP::Mesh::Vertex, "Normal");
    siloWriter->registerVector(suckItVec, meshAdapter, AMP::Mesh::Vertex, "Suction");
    siloWriter->registerVector(lickItVec, meshAdapter, AMP::Mesh::Vertex, "Tangent");
    siloWriter->registerVector(contactShiftVec, meshAdapter, AMP::Mesh::Vertex, "Shift");
    char outFileName[256];
    sprintf(outFileName, "TITI_%d", 0);
    siloWriter->writeFile(outFileName, 0);
  }
#endif

//  bool skipDisplaceMesh = true;
//  contactOperator->updateActiveSet(nullVec, skipDisplaceMesh);

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
  AMP_ASSERT( matLocalSize == static_cast<int>(dispDofManager->numLocalDOF()) );
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
//  linearSolver->setZeroInitialGuess(true);
  linearSolver->setInitialGuess(columnSolVec);

  std::vector<AMP::Mesh::MeshElementID> slaveNodesGlobalIDs;
  selectNodes(slaveMeshAdapter, slaveNodesGlobalIDs);
  printNodesValues(slaveMeshAdapter, slaveNodesGlobalIDs, suckItVec, lickItVec);

int TOTO_count = 0;
size_t const maxLoadingIterations = input_db->getIntegerWithDefault("maxLoadingIterations", 5);
for (size_t loadingIteration = 0; loadingIteration < maxLoadingIterations; ++loadingIteration) {
  double scalingFactor = static_cast<double>(loadingIteration+1) / static_cast<double>(maxLoadingIterations);
  // quadratic scheme
//  scalingFactor = scalingFactor * scalingFactor;
  if (!rank) { std::cout<<"LOADING STEP "<<loadingIteration+1<<"/"<<maxLoadingIterations<<"  (factor="<<scalingFactor<<")\n"; }

  size_t const maxActiveSetIterations = input_db->getIntegerWithDefault("maxActiveSetIterations", 5);
  for (size_t activeSetIteration = 0; activeSetIteration < maxActiveSetIterations; ++activeSetIteration) {
  if (!rank) { std::cout<<"ACTIVE SET ITERATION #"<<activeSetIteration+1<<"\n"; }
++TOTO_count;

  columnSolVec->zero();
  columnRhsVec->zero();

  double masterLoadParameter = input_db->getDouble("MasterLoadParameter");
  double masterLoadCutoff = input_db->getDouble("MasterLoadCutoff");
  getConcentratedLoadAtNodes(masterLoadParameter, masterLoadCutoff, masterMeshAdapter, columnRhsVec, dispDofManager);
  columnRhsVec->scale(scalingFactor);

/*
if (masterMeshAdapter.get() != NULL) {
AMP::Mesh::MeshIterator masterBoundaryIterator = masterMeshAdapter->getBoundaryIDIterator(AMP::Mesh::Vertex, 4, 0);
AMP::Mesh::MeshIterator masterBoundaryIterator_begin = masterBoundaryIterator.begin(),
    masterBoundaryIterator_end = masterBoundaryIterator.end();
std::vector<double> vertexCoordinates;

size_t count = 0;
for (masterBoundaryIterator = masterBoundaryIterator_begin; masterBoundaryIterator != masterBoundaryIterator_end; ++masterBoundaryIterator) {
  vertexCoordinates = masterBoundaryIterator->coord();
  AMP_ASSERT( vertexCoordinates.size() == 3 );
  fout<<std::setprecision(15);
  double masterLoadParameter = input_db->getDouble("SlaveLoadParameter");
  if (vertexCoordinates[1] > input_db->getDouble("SlaveLoadCutoff")) {
    ++count;
    fout<<vertexCoordinates[0]<<"  "<<vertexCoordinates[1]<<"  "<<vertexCoordinates[2]<<std::endl;
    std::vector<size_t> DOFsIndices;
    dispDofManager->getDOFs(masterBoundaryIterator->globalID(), DOFsIndices);
    AMP_ASSERT( DOFsIndices.size() == 3 );
    columnRhsVec->setLocalValueByGlobalID(DOFsIndices[1], masterLoadParameter);
  } // end if
} // end for
size_t count_total = masterMeshAdapter->getComm().sumReduce(count);
if (!masterMeshAdapter->getComm().getRank()) {
std::cout<<"count_total="<<count_total<<"\n";
}
columnRhsVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
} // end if
*/


  // apply dirichlet rhs correction
  if (masterBVPOperator.get() != NULL) {
    masterBVPOperator->modifyRHSvector(columnRhsVec);
  } // end if
  if (slaveBVPOperator.get() != NULL) {
    slaveBVPOperator->modifyRHSvector(columnRhsVec);
  } // end if

  // get d
  contactOperator->addShiftToSlave(columnSolVec);

  // compute - Kd
  AMP::LinearAlgebra::Vector::shared_ptr rhsCorrectionVec = createVector(dispDofManager, columnVar, split);
  columnOperator->apply(nullVec, columnSolVec, rhsCorrectionVec, -1.0, 0.0);
  columnOperator->append(contactOperator);

  // f = f - Kd
  columnRhsVec->add(columnRhsVec, rhsCorrectionVec);

  // f^m = f^m + C^T f^s
  // f^s = 0
  contactOperator->addSlaveToMaster(columnRhsVec);
  contactOperator->setSlaveToZero(columnRhsVec);

  // u_s = C u_m
  contactOperator->copyMasterToSlave(columnSolVec);

  linearSolver->solve(columnRhsVec, columnSolVec);

  // u^s = C u^m + d
  contactOperator->copyMasterToSlave(columnSolVec);
  contactOperator->addShiftToSlave(columnSolVec);

  computeStressTensor(masterMeshAdapter, columnSolVec, 
      sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy,
      sigma_eff, masterMechanicsMaterialModel,
      referenceTemperature, thermalExpansionCoefficient, tempVec);
  computeStressTensor(slaveMeshAdapter, columnSolVec, 
      sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy,
      sigma_eff, slaveMechanicsMaterialModel,
      referenceTemperature, thermalExpansionCoefficient, tempVec);

  activeSetVec->setToScalar(-1.0);

  std::vector<AMP::Mesh::MeshElementID> const * pointerToActiveSet;
  contactOperator->getActiveSet(pointerToActiveSet);
  size_t const sizeOfActiveSetBeforeUpdate = pointerToActiveSet->size();

  std::vector<size_t> activeSetTempDOFsIndicesBeforeUpdate;
  tempDofManager->getDOFs(*pointerToActiveSet, activeSetTempDOFsIndicesBeforeUpdate);
  AMP_ASSERT( activeSetTempDOFsIndicesBeforeUpdate.size() == sizeOfActiveSetBeforeUpdate );
  std::vector<double> valuesForActiveSet(sizeOfActiveSetBeforeUpdate, 2.0); 
  activeSetVec->setLocalValuesByGlobalID(sizeOfActiveSetBeforeUpdate, &(activeSetTempDOFsIndicesBeforeUpdate[0]), &(valuesForActiveSet[0]));

  std::vector<size_t> activeSetDispDOFsIndicesBeforeUpdate;
  dispDofManager->getDOFs(*pointerToActiveSet, activeSetDispDOFsIndicesBeforeUpdate);
  AMP_ASSERT( activeSetDispDOFsIndicesBeforeUpdate.size() == 3*sizeOfActiveSetBeforeUpdate );

  {//TODO: this is tmp only
  normalVectorVec->zero();
  std::vector<double> const * slaveVerticesNormalVector;
  std::vector<double> const * slaveVerticesSurfaceTraction;
  contactOperator->getSlaveVerticesNormalVectorAndSurfaceTraction(slaveVerticesNormalVector, slaveVerticesSurfaceTraction);
  AMP_ASSERT( slaveVerticesNormalVector->size() == 3*sizeOfActiveSetBeforeUpdate );
//  AMP_ASSERT( slaveVerticesSurfaceTraction->size() == 3*sizeOfActiveSetBeforeUpdate );
  normalVectorVec->setLocalValuesByGlobalID(3*sizeOfActiveSetBeforeUpdate, &(activeSetDispDOFsIndicesBeforeUpdate[0]), &((*slaveVerticesNormalVector)[0]));
  }

#ifdef USE_EXT_SILO
  {
  meshAdapter->displaceMesh(columnSolVec);
  char outFileName[256];
  sprintf(outFileName, "TITI_%d", 0);
//  siloWriter->writeFile(outFileName, activeSetIteration+1);
  siloWriter->writeFile(outFileName, TOTO_count);
  columnSolVec->scale(-1.0);
  meshAdapter->displaceMesh(columnSolVec);
  columnSolVec->scale(-1.0);
  }
#endif

//  meshAdapter->displaceMesh(columnSolVec);
  size_t nChangesInActiveSet = contactOperator->updateActiveSet(columnSolVec);
  if (!rank) { std::cout<<nChangesInActiveSet<<" CHANGES IN ACTIVE SET\n"; }

  suckItVec->zero();
  lickItVec->zero();
  surfaceTractionVec->zero();
  normalVectorVec->zero();
  size_t const sizeOfActiveSetAfterUpdate = pointerToActiveSet->size();
  std::vector<double> const * slaveVerticesNormalVector;
  std::vector<double> const * slaveVerticesSurfaceTraction;
  contactOperator->getSlaveVerticesNormalVectorAndSurfaceTraction(slaveVerticesNormalVector, slaveVerticesSurfaceTraction);
//  AMP_ASSERT( slaveVerticesNormalVector->size() == 3*sizeOfActiveSetBeforeUpdate );
  AMP_ASSERT( slaveVerticesNormalVector->size() == 3*sizeOfActiveSetAfterUpdate );
  AMP_ASSERT( slaveVerticesSurfaceTraction->size() == 3*sizeOfActiveSetBeforeUpdate );
  surfaceTractionVec->setLocalValuesByGlobalID(3*sizeOfActiveSetBeforeUpdate, &(activeSetDispDOFsIndicesBeforeUpdate[0]), &((*slaveVerticesSurfaceTraction)[0]));
  normalVectorVec->setLocalValuesByGlobalID(3*sizeOfActiveSetBeforeUpdate, &(activeSetDispDOFsIndicesBeforeUpdate[0]), &((*slaveVerticesNormalVector)[0]));
  std::vector<double> surfaceTractionDOTnormalVector(sizeOfActiveSetBeforeUpdate);
  std::vector<double> surfaceTractionTangentComponent(sizeOfActiveSetBeforeUpdate);
  for (size_t kk = 0; kk < sizeOfActiveSetBeforeUpdate;++kk) {
    surfaceTractionDOTnormalVector[kk] = compute_scalar_product(&((*slaveVerticesSurfaceTraction)[3*kk]), &((*slaveVerticesNormalVector)[3*kk]));
    surfaceTractionTangentComponent[kk] = std::sqrt(std::pow(compute_vector_norm(&((*slaveVerticesSurfaceTraction)[3*kk])), 2) - std::pow(surfaceTractionDOTnormalVector[kk], 2));
  } // end for kk
  suckItVec->setLocalValuesByGlobalID(sizeOfActiveSetBeforeUpdate, &(activeSetTempDOFsIndicesBeforeUpdate[0]), &(surfaceTractionDOTnormalVector[0]));
  lickItVec->setLocalValuesByGlobalID(sizeOfActiveSetBeforeUpdate, &(activeSetTempDOFsIndicesBeforeUpdate[0]), &(surfaceTractionTangentComponent[0]));

  printNodesValues(slaveMeshAdapter, slaveNodesGlobalIDs, suckItVec, lickItVec);

#ifdef USE_EXT_SILO
  {
  meshAdapter->displaceMesh(columnSolVec);
  char outFileName[256];
  sprintf(outFileName, "TITI_%d", 0);
//  siloWriter->writeFile(outFileName, activeSetIteration+1);
  siloWriter->writeFile(outFileName, TOTO_count);
  columnSolVec->scale(-1.0);
  meshAdapter->displaceMesh(columnSolVec);
  columnSolVec->scale(-1.0);
  }
#endif

  if (nChangesInActiveSet == 0) { break; }
  AMP_ASSERT( activeSetIteration != maxActiveSetIterations - 1 );
  } // end for

} //end for

  meshAdapter->displaceMesh(columnSolVec);

if (masterMeshAdapter.get() != NULL) {
std::fstream masterFout;
masterFout.open("master_pellet_displaced_mesh", std::fstream::out);
double point_of_view[3] = { 1.0, 1.0, 1.0 };
drawFacesOnBoundaryID(masterMeshAdapter, 1, masterFout, point_of_view, "");
drawFacesOnBoundaryID(masterMeshAdapter, 4, masterFout, point_of_view, "");
masterFout.close();
} // end if
if (slaveMeshAdapter.get() != NULL) {
std::fstream slaveFout;
slaveFout.open("slave_pellet_displaced_mesh", std::fstream::out);
double point_of_view[3] = { 1.0, 1.0, 1.0 };
drawFacesOnBoundaryID(slaveMeshAdapter, 1, slaveFout, point_of_view, "dashed");
drawFacesOnBoundaryID(slaveMeshAdapter, 4, slaveFout, point_of_view, "dashed");
//drawVerticesOnBoundaryID(slaveMeshAdapter, 2, slaveFout, point_of_view, "red");
slaveFout.close();
} // end if

#ifdef USE_EXT_SILO
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
  AMP::UnitTest ut;

  std::vector<std::string> exeNames; 
  exeNames.push_back("testNodeToFaceContactOperator-2");

  try {
    for (size_t i = 0; i < exeNames.size(); ++i) { 
      myTest(&ut, exeNames[i]); 
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



