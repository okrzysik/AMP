
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
#include "operators/TrilinosMatrixShellOperator.h"
#include "operators/petsc/PetscMatrixShellOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/mechanics/MechanicsModelParameters.h"
#include "operators/mechanics/MechanicsMaterialModel.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/ConstructLinearMechanicsRHSVector.h"
#include "operators/contact/NodeToFaceContactOperator.h"
#include "operators/mechanics/IsotropicElasticModel.h"

#include "solvers/ColumnSolver.h"
#include "solvers/ConstraintsEliminationSolver.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

#include "utils/ReadTestMesh.h"

#include <set>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include "ampmesh/latex_visualization_tools.h"
#include "ampmesh/euclidean_geometry_tools.h"

#include "testNodeToFaceContactOperator.h"

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
  boost::shared_ptr<AMP::Solver::ColumnSolverParameters> columnPreconditionerParams(new AMP::Solver::ColumnSolverParameters(columnPreconditioner_db));
  columnPreconditionerParams->d_pOperator = columnOperator;
  boost::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

  // Get the mechanics material model for the contact operator
  boost::shared_ptr<AMP::Database> model_db = input_db->getDatabase("MechanicsMaterialModel");
  boost::shared_ptr<AMP::Operator::MechanicsModelParameters> mechanicsMaterialModelParams(new AMP::Operator::MechanicsModelParameters(model_db));
  boost::shared_ptr<AMP::Operator::MechanicsMaterialModel> masterMechanicsMaterialModel(new AMP::Operator::IsotropicElasticModel(mechanicsMaterialModelParams));

  // Build the contact operators
  boost::shared_ptr<AMP::Operator::NodeToFaceContactOperator> bottomPelletTopPelletContactOperator;
  boost::shared_ptr<AMP::Operator::NodeToFaceContactOperator> bottomPelletCladContactOperator;
  boost::shared_ptr<AMP::Operator::NodeToFaceContactOperator> topPelletCladContactOperator;

  boost::shared_ptr<AMP::Database> bottomPelletTopPelletContact_db = input_db->getDatabase("BottomPelletTopPelletContactOperator");
  boost::shared_ptr<AMP::Operator::ContactOperatorParameters> bottomPelletTopPelletContactOperatorParams(new AMP::Operator::ContactOperatorParameters(bottomPelletTopPelletContact_db));
  bottomPelletTopPelletContactOperatorParams->d_DOFsPerNode = dofsPerNode;
  bottomPelletTopPelletContactOperatorParams->d_DOFManager = dispDofManager;
  bottomPelletTopPelletContactOperatorParams->d_GlobalComm = globalComm;
  bottomPelletTopPelletContactOperatorParams->d_Mesh = meshAdapter;
  bottomPelletTopPelletContactOperatorParams->d_MasterMechanicsMaterialModel = masterMechanicsMaterialModel;
  bottomPelletTopPelletContactOperatorParams->reset();
  bottomPelletTopPelletContactOperator = boost::shared_ptr<AMP::Operator::NodeToFaceContactOperator>(new AMP::Operator::NodeToFaceContactOperator(bottomPelletTopPelletContactOperatorParams));

  boost::shared_ptr<AMP::Database> bottomPelletCladContact_db = input_db->getDatabase("BottomPelletCladContactOperator");
  boost::shared_ptr<AMP::Operator::ContactOperatorParameters> bottomPelletCladContactOperatorParams(new AMP::Operator::ContactOperatorParameters(bottomPelletCladContact_db));
  bottomPelletCladContactOperatorParams->d_DOFsPerNode = dofsPerNode;
  bottomPelletCladContactOperatorParams->d_DOFManager = dispDofManager;
  bottomPelletCladContactOperatorParams->d_GlobalComm = globalComm;
  bottomPelletCladContactOperatorParams->d_Mesh = meshAdapter;
  bottomPelletCladContactOperatorParams->d_MasterMechanicsMaterialModel = masterMechanicsMaterialModel;
  bottomPelletCladContactOperatorParams->reset();
  bottomPelletCladContactOperator = boost::shared_ptr<AMP::Operator::NodeToFaceContactOperator>(new AMP::Operator::NodeToFaceContactOperator(bottomPelletCladContactOperatorParams));

  boost::shared_ptr<AMP::Database> topPelletCladContact_db = input_db->getDatabase("TopPelletCladContactOperator");
  boost::shared_ptr<AMP::Operator::ContactOperatorParameters> topPelletCladContactOperatorParams(new AMP::Operator::ContactOperatorParameters(topPelletCladContact_db));
  topPelletCladContactOperatorParams->d_DOFsPerNode = dofsPerNode;
  topPelletCladContactOperatorParams->d_DOFManager = dispDofManager;
  topPelletCladContactOperatorParams->d_GlobalComm = globalComm;
  topPelletCladContactOperatorParams->d_Mesh = meshAdapter;
  topPelletCladContactOperatorParams->d_MasterMechanicsMaterialModel = masterMechanicsMaterialModel;
  topPelletCladContactOperatorParams->reset();
  topPelletCladContactOperator = boost::shared_ptr<AMP::Operator::NodeToFaceContactOperator>(new AMP::Operator::NodeToFaceContactOperator(topPelletCladContactOperatorParams));

  bottomPelletTopPelletContactOperator->initialize();
  bottomPelletCladContactOperator->initialize();
  topPelletCladContactOperator->initialize();
  
  bool useML = input_db->getBoolWithDefault("useML", false);
  // Build the BVP operators
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> bottomPelletBVPOperator;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> topPelletBVPOperator;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> cladBVPOperator;


  AMP::Mesh::MeshID bottomPelletMeshID = bottomPelletTopPelletContactOperator->getSlaveMeshID();
  AMP_ASSERT(bottomPelletMeshID == bottomPelletCladContactOperator->getSlaveMeshID());
  AMP::Mesh::Mesh::shared_ptr bottomPelletMeshAdapter = meshAdapter->Subset(bottomPelletMeshID);
  if (bottomPelletMeshAdapter.get() != NULL) {
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> bottomPelletElementPhysicsModel;
    bottomPelletBVPOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(bottomPelletMeshAdapter,
                                                                                                                                           "BottomPelletBVPOperator",
                                                                                                                                           input_db,
                                                                                                                                           bottomPelletElementPhysicsModel));
    columnOperator->append(bottomPelletBVPOperator);

    boost::shared_ptr<AMP::Database> bottomPelletSolver_db = columnPreconditioner_db->getDatabase("DummySolver"); 
    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> bottomPelletSolverParams(new AMP::Solver::PetscKrylovSolverParameters(bottomPelletSolver_db));
    bottomPelletSolverParams->d_pOperator = bottomPelletBVPOperator;
    bottomPelletSolverParams->d_comm = bottomPelletMeshAdapter->getComm();
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> bottomPelletSolver(new AMP::Solver::PetscKrylovSolver(bottomPelletSolverParams));
    columnPreconditioner->append(bottomPelletSolver);
  } // end if

  AMP::Mesh::MeshID topPelletMeshID = bottomPelletTopPelletContactOperator->getMasterMeshID();
  AMP_ASSERT(topPelletMeshID == topPelletCladContactOperator->getSlaveMeshID());
  AMP::Mesh::Mesh::shared_ptr topPelletMeshAdapter = meshAdapter->Subset(topPelletMeshID);
  if (topPelletMeshAdapter.get() != NULL) {
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> topPelletElementPhysicsModel;
    topPelletBVPOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(topPelletMeshAdapter,
                                                                                                                                        "TopPelletBVPOperator",
                                                                                                                                        input_db,
                                                                                                                                        topPelletElementPhysicsModel));
    columnOperator->append(topPelletBVPOperator);

    boost::shared_ptr<AMP::Database> topPelletSolver_db = columnPreconditioner_db->getDatabase("DummySolver"); 
    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> topPelletSolverParams(new AMP::Solver::PetscKrylovSolverParameters(topPelletSolver_db));
    topPelletSolverParams->d_pOperator = topPelletBVPOperator;
    topPelletSolverParams->d_comm = topPelletMeshAdapter->getComm();
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> topPelletSolver(new AMP::Solver::PetscKrylovSolver(topPelletSolverParams));
    columnPreconditioner->append(topPelletSolver);
  } // end if

  AMP::Mesh::MeshID cladMeshID = bottomPelletCladContactOperator->getMasterMeshID();
  AMP_ASSERT(cladMeshID == topPelletCladContactOperator->getMasterMeshID());
  AMP::Mesh::Mesh::shared_ptr cladMeshAdapter = meshAdapter->Subset(cladMeshID);
  if (cladMeshAdapter.get() != NULL) {
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> cladElementPhysicsModel;
    cladBVPOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(cladMeshAdapter,
                                                                                                                                   "CladBVPOperator",
                                                                                                                                   input_db,
                                                                                                                                   cladElementPhysicsModel));
    columnOperator->append(cladBVPOperator);

    if (useML) {
      boost::shared_ptr<AMP::Database> cladSolver_db = columnPreconditioner_db->getDatabase("MLSolver"); 
      boost::shared_ptr<AMP::Solver::SolverStrategyParameters> cladSolverParams(new AMP::Solver::SolverStrategyParameters(cladSolver_db));
      cladSolverParams->d_pOperator = cladBVPOperator;
      boost::shared_ptr<AMP::Solver::TrilinosMLSolver> cladSolver(new AMP::Solver::TrilinosMLSolver(cladSolverParams));
      columnPreconditioner->append(cladSolver);
    } else {
      boost::shared_ptr<AMP::Database> cladSolver_db = columnPreconditioner_db->getDatabase("DummySolver"); 
      boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> cladSolverParams(new AMP::Solver::PetscKrylovSolverParameters(cladSolver_db));
      cladSolverParams->d_pOperator = cladBVPOperator;
      cladSolverParams->d_comm = cladMeshAdapter->getComm();
      boost::shared_ptr<AMP::Solver::PetscKrylovSolver> cladSolver(new AMP::Solver::PetscKrylovSolver(cladSolverParams));
      columnPreconditioner->append(cladSolver);
   } // end if

  } // end if


{
  boost::shared_ptr<AMP::Database> bottomPelletTopPelletContactPreconditioner_db = columnPreconditioner_db->getDatabase("ContactPreconditioner"); 
  boost::shared_ptr<AMP::Solver::ConstraintsEliminationSolverParameters> bottomPelletTopPelletContactPreconditionerParams(new AMP::Solver::ConstraintsEliminationSolverParameters(bottomPelletTopPelletContactPreconditioner_db));
  bottomPelletTopPelletContactPreconditionerParams->d_pOperator = bottomPelletTopPelletContactOperator;
  boost::shared_ptr<AMP::Solver::ConstraintsEliminationSolver> bottomPelletTopPelletContactPreconditioner(new AMP::Solver::ConstraintsEliminationSolver(bottomPelletTopPelletContactPreconditionerParams));
  columnPreconditioner->append(bottomPelletTopPelletContactPreconditioner);
}
{
  boost::shared_ptr<AMP::Database> bottomPelletCladContactPreconditioner_db = columnPreconditioner_db->getDatabase("ContactPreconditioner"); 
  boost::shared_ptr<AMP::Solver::ConstraintsEliminationSolverParameters> bottomPelletCladContactPreconditionerParams(new AMP::Solver::ConstraintsEliminationSolverParameters(bottomPelletCladContactPreconditioner_db));
  bottomPelletCladContactPreconditionerParams->d_pOperator = bottomPelletCladContactOperator;
  boost::shared_ptr<AMP::Solver::ConstraintsEliminationSolver> bottomPelletCladContactPreconditioner(new AMP::Solver::ConstraintsEliminationSolver(bottomPelletCladContactPreconditionerParams));
  columnPreconditioner->append(bottomPelletCladContactPreconditioner);
}
{
  boost::shared_ptr<AMP::Database> topPelletCladContactPreconditioner_db = columnPreconditioner_db->getDatabase("ContactPreconditioner"); 
  boost::shared_ptr<AMP::Solver::ConstraintsEliminationSolverParameters> topPelletCladContactPreconditionerParams(new AMP::Solver::ConstraintsEliminationSolverParameters(topPelletCladContactPreconditioner_db));
  topPelletCladContactPreconditionerParams->d_pOperator = topPelletCladContactOperator;
  boost::shared_ptr<AMP::Solver::ConstraintsEliminationSolver> topPelletCladContactPreconditioner(new AMP::Solver::ConstraintsEliminationSolver(topPelletCladContactPreconditionerParams));
  columnPreconditioner->append(topPelletCladContactPreconditioner);
}

  // Items for computing the RHS correction due to thermal expansion
  boost::shared_ptr<AMP::Database> temperatureRhs_db = input_db->getDatabase("TemperatureRHSVectorCorrection");
  AMP::LinearAlgebra::Variable::shared_ptr tempVar(new AMP::LinearAlgebra::Variable("temperature"));
  AMP::LinearAlgebra::Variable::shared_ptr dispVar = columnOperator->getOutputVariable();
  AMP::Discretization::DOFManager::shared_ptr tempDofManager = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth, 1 , split);
  AMP::LinearAlgebra::Vector::shared_ptr tempVec = AMP::LinearAlgebra::createVector(tempDofManager, tempVar, split);
  AMP::LinearAlgebra::Vector::shared_ptr refTempVec = tempVec->cloneVector();

  AMP::LinearAlgebra::Vector::shared_ptr sigma_xx = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_xx")), split);
  AMP::LinearAlgebra::Vector::shared_ptr sigma_yy = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_yy")), split);
  AMP::LinearAlgebra::Vector::shared_ptr sigma_zz = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_zz")), split);
  AMP::LinearAlgebra::Vector::shared_ptr sigma_yz = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_yz")), split);
  AMP::LinearAlgebra::Vector::shared_ptr sigma_xz = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_xz")), split);
  AMP::LinearAlgebra::Vector::shared_ptr sigma_xy = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_xy")), split);
  AMP::LinearAlgebra::Vector::shared_ptr sigma_eff = AMP::LinearAlgebra::createVector(tempDofManager, AMP::LinearAlgebra::Variable::shared_ptr(new AMP::LinearAlgebra::Variable("sigma_eff")), split);

  AMP::Mesh::MeshIterator topPelletMeshIterator = topPelletMeshAdapter->getIterator(AMP::Mesh::Vertex);
  AMP::Mesh::MeshIterator bottomPelletMeshIterator = bottomPelletMeshAdapter->getIterator(AMP::Mesh::Vertex);
  AMP::Mesh::MeshIterator meshIterator = AMP::Mesh::Mesh::getIterator(AMP::Mesh::Union, topPelletMeshIterator, bottomPelletMeshIterator),
      meshIterator_begin = meshIterator.begin(),
      meshIterator_end = meshIterator.end();
  std::vector<double> vertexCoord;
  std::vector<size_t> DOFsIndices;
  double temperatureOuterRadius = input_db->getDouble("TemperatureOuterRadius"); 
  double heatGenerationRate = input_db->getDouble("HeatGenerationRate");
  double outerRadius = input_db->getDouble("OuterRadius");
  double outerRadiusSquared = outerRadius * outerRadius;
  double thermalConductivity = input_db->getDouble("ThermalConductivity");
  double temperatureCenterLine = temperatureOuterRadius + heatGenerationRate * outerRadiusSquared / (4.0 * thermalConductivity);
  double referenceTemperature = input_db->getDouble("ReferenceTemperature");
  if (!rank) { std::cout<<"temperatureCenterLine="<<temperatureCenterLine<<"\n"; }
  refTempVec->setToScalar(referenceTemperature);
  tempVec->setToScalar(0.0);

  for (meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator) {
    vertexCoord = meshIterator->coord();
//rotate_points(2, M_PI / -2.0, 1, &(vertexCoord[0]));
    double radiusSquared = vertexCoord[0]*vertexCoord[0] + vertexCoord[1]*vertexCoord[1];
    double temperature = temperatureCenterLine - heatGenerationRate * radiusSquared / (4.0 * thermalConductivity);
    tempDofManager->getDOFs(meshIterator->globalID(), DOFsIndices);
    AMP_ASSERT(DOFsIndices.size() == 1);
    tempVec->setLocalValuesByGlobalID(1, &(DOFsIndices[0]), &temperature);
  } // end for

  AMP::LinearAlgebra::VS_Mesh cladVectorSelector(cladMeshAdapter);
  AMP::LinearAlgebra::Vector::shared_ptr cladTempVec = tempVec->select(cladVectorSelector, tempVar->getName());
  cladTempVec->setToScalar(referenceTemperature);

  boost::shared_ptr<AMP::Database> tmp_db = temperatureRhs_db->getDatabase("RhsMaterialModel");
  double thermalExpansionCoefficient = tmp_db->getDouble("THERMAL_EXPANSION_COEFFICIENT");
  bottomPelletTopPelletContactOperator->uglyHack(tempVec, tempDofManager, thermalExpansionCoefficient, referenceTemperature);
  bottomPelletCladContactOperator->uglyHack(tempVec, tempDofManager, thermalExpansionCoefficient, referenceTemperature);
  topPelletCladContactOperator->uglyHack(tempVec, tempDofManager, thermalExpansionCoefficient, referenceTemperature);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = AMP::LinearAlgebra::createVector(dispDofManager, columnVar, split);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = AMP::LinearAlgebra::createVector(dispDofManager, columnVar, split);
  columnSolVec->zero();
  columnRhsVec->zero();

  AMP::LinearAlgebra::Vector::shared_ptr activeSetVec = sigma_eff->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr suckItVec = sigma_eff->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr surfaceTractionVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr normalVectorVec = columnSolVec->cloneVector();

  computeStressTensor(meshAdapter, columnSolVec, 
      sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy,
      sigma_eff, 1.0e6, 0.3,
      referenceTemperature, thermalExpansionCoefficient, tempVec);

  bool skipDisplaceMesh = true;
  bottomPelletTopPelletContactOperator->updateActiveSet(nullVec, skipDisplaceMesh);
  bottomPelletCladContactOperator->updateActiveSet(nullVec, skipDisplaceMesh);
  topPelletCladContactOperator->updateActiveSet(nullVec, skipDisplaceMesh);

  
  AMP::LinearAlgebra::Vector::shared_ptr contactShiftVec = createVector(dispDofManager, columnVar, split);
  contactShiftVec->zero();

  AMP::LinearAlgebra::Vector::shared_ptr oldSolVec = columnSolVec->cloneVector();
  oldSolVec->zero();

#ifdef USE_EXT_SILO
  {
    siloWriter->registerVector(columnSolVec, meshAdapter, AMP::Mesh::Vertex, "Solution");
    siloWriter->registerVector(tempVec, meshAdapter, AMP::Mesh::Vertex, "Temperature");
    siloWriter->registerVector(sigma_eff, meshAdapter, AMP::Mesh::Vertex, "vonMises");
    siloWriter->registerVector(sigma_xx, meshAdapter, AMP::Mesh::Vertex, "sigma_xx");
    siloWriter->registerVector(sigma_yy, meshAdapter, AMP::Mesh::Vertex, "sigma_yy");
    siloWriter->registerVector(sigma_zz, meshAdapter, AMP::Mesh::Vertex, "sigma_zz");
    siloWriter->registerVector(sigma_yz, meshAdapter, AMP::Mesh::Vertex, "sigma_yz");
    siloWriter->registerVector(sigma_xz, meshAdapter, AMP::Mesh::Vertex, "sigma_xz");
    siloWriter->registerVector(sigma_xy, meshAdapter, AMP::Mesh::Vertex, "sigma_xy");
    siloWriter->registerVector(activeSetVec, meshAdapter, AMP::Mesh::Vertex, "Contact");
    siloWriter->registerVector(oldSolVec, meshAdapter, AMP::Mesh::Vertex, "Error");
    siloWriter->registerVector(surfaceTractionVec, meshAdapter, AMP::Mesh::Vertex, "Traction");
    siloWriter->registerVector(normalVectorVec, meshAdapter, AMP::Mesh::Vertex, "Normal");
    siloWriter->registerVector(suckItVec, meshAdapter, AMP::Mesh::Vertex, "Suction");
    siloWriter->registerVector(contactShiftVec, meshAdapter, AMP::Mesh::Vertex, "Shift");
    char outFileName[256];
    sprintf(outFileName, "TATA_%d", 0);
    siloWriter->writeFile(outFileName, 0);
  }
#endif
  oldSolVec->copyVector(columnSolVec);

  columnSolVec->zero();
  columnOperator->append(bottomPelletTopPelletContactOperator);
  columnOperator->append(bottomPelletCladContactOperator);
  columnOperator->append(topPelletCladContactOperator);


  // Build a matrix shell operator to use the column operator with the petsc krylov solvers
  boost::shared_ptr<AMP::Database> matrixShellDatabase = input_db->getDatabase("MatrixShellOperator");
  boost::shared_ptr<AMP::Operator::OperatorParameters> matrixShellParams(new AMP::Operator::OperatorParameters(matrixShellDatabase));
  boost::shared_ptr<AMP::Operator::PetscMatrixShellOperator> matrixShellOperator(new AMP::Operator::PetscMatrixShellOperator(matrixShellParams));

  int numBottomPelletLocalNodes = 0;
  int numTopPelletLocalNodes = 0;
  int numCladLocalNodes = 0;
  if (bottomPelletMeshAdapter.get() != NULL) { numBottomPelletLocalNodes = bottomPelletMeshAdapter->numLocalElements(AMP::Mesh::Vertex); }
  if (topPelletMeshAdapter.get() != NULL) { numTopPelletLocalNodes = topPelletMeshAdapter->numLocalElements(AMP::Mesh::Vertex); }
  if (cladMeshAdapter.get() != NULL) { numCladLocalNodes = cladMeshAdapter->numLocalElements(AMP::Mesh::Vertex); }
  int matLocalSize = dofsPerNode * (numBottomPelletLocalNodes + numTopPelletLocalNodes + numCladLocalNodes);
  AMP_ASSERT( matLocalSize == static_cast<int>(dispDofManager->numLocalDOF()) );
  matrixShellOperator->setComm(globalComm);
  matrixShellOperator->setMatLocalRowSize(matLocalSize);
  matrixShellOperator->setMatLocalColumnSize(matLocalSize);
  matrixShellOperator->setOperator(columnOperator); 

  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = matrixShellOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = columnPreconditioner;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));
//  linearSolver->setZeroInitialGuess(true);
  linearSolver->setInitialGuess(columnSolVec);

  AMP::LinearAlgebra::Vector::shared_ptr fullThermalLoadingTempMinusRefTempVec = tempVec->cloneVector();
  fullThermalLoadingTempMinusRefTempVec->subtract(tempVec, refTempVec);

  size_t maxActiveSetIterations = input_db->getIntegerWithDefault("maxActiveSetIterations", 5);
  size_t maxThermalLoadingIterations = input_db->getIntegerWithDefault("maxThermalLoadingIterations", 5);
  std::vector<double> scalingFactors;
  std::vector<int> maxIterations;
  bool customLoading = input_db->getBoolWithDefault("customLoading", false);
  if (customLoading) {
    scalingFactors = input_db->getDoubleArray("scalingFactors");
    maxIterations = input_db->getIntegerArray("maxIterations");
    AMP_ASSERT(scalingFactors.size() == maxIterations.size());
    maxThermalLoadingIterations = scalingFactors.size();
  } // end if

int TOTO_count = 0;
for (size_t thermalLoadingIteration = 0; thermalLoadingIteration < maxThermalLoadingIterations; ++thermalLoadingIteration) {
  double scalingFactor = static_cast<double>(thermalLoadingIteration+1) / static_cast<double>(maxThermalLoadingIterations);
  if (customLoading) {
    scalingFactor = scalingFactors[thermalLoadingIteration];
    maxActiveSetIterations = maxIterations[thermalLoadingIteration];
  } // end if

  if (!rank) { std::cout<<"THERMAL LOADING "<<thermalLoadingIteration+1<<"/"<<maxThermalLoadingIterations<<"  ("<<scalingFactor<<")\n"; }
  tempVec->axpy(scalingFactor, fullThermalLoadingTempMinusRefTempVec, refTempVec);

  for (size_t activeSetIteration = 0; activeSetIteration < maxActiveSetIterations; ++activeSetIteration) {
    if (!rank) { std::cout<<"ACTIVE SET ITERATION #"<<activeSetIteration+1<<std::endl; }
++ TOTO_count;

    columnSolVec->zero();
    columnRhsVec->zero();

    // compute thermal load f
    computeTemperatureRhsVector(meshAdapter, temperatureRhs_db, tempVar, dispVar, tempVec, refTempVec, columnRhsVec);

    // apply dirichlet rhs correction on f
    if (bottomPelletBVPOperator.get() != NULL) { bottomPelletBVPOperator->modifyRHSvector(columnRhsVec); }
    if (topPelletBVPOperator.get() != NULL) { topPelletBVPOperator->modifyRHSvector(columnRhsVec); }
    if (cladBVPOperator.get() != NULL) { cladBVPOperator->modifyRHSvector(columnRhsVec); }

    // get d
    contactShiftVec->zero();
    bottomPelletTopPelletContactOperator->addShiftToSlave(contactShiftVec);
    bottomPelletCladContactOperator->addShiftToSlave(contactShiftVec);
    topPelletCladContactOperator->addShiftToSlave(contactShiftVec);

    // compute - Kd
    AMP::LinearAlgebra::Vector::shared_ptr rhsCorrectionVec = createVector(dispDofManager, columnVar, split);
    rhsCorrectionVec->zero();
    bottomPelletBVPOperator->apply(nullVec, contactShiftVec, rhsCorrectionVec, -1.0, 0.0);
    topPelletBVPOperator->apply(nullVec, contactShiftVec, rhsCorrectionVec, -1.0, 0.0);
    cladBVPOperator->apply(nullVec, contactShiftVec, rhsCorrectionVec, -1.0, 0.0);

    // f = f - Kd
    columnRhsVec->add(columnRhsVec, rhsCorrectionVec);

    // f^m = f^m + C^T f^s
    // f^s = 0
    bottomPelletTopPelletContactOperator->addSlaveToMaster(columnRhsVec);
    bottomPelletCladContactOperator->addSlaveToMaster(columnRhsVec);
    topPelletCladContactOperator->addSlaveToMaster(columnRhsVec);

    bottomPelletTopPelletContactOperator->setSlaveToZero(columnRhsVec);
    bottomPelletCladContactOperator->setSlaveToZero(columnRhsVec);
    topPelletCladContactOperator->setSlaveToZero(columnRhsVec);

    // u_s = C u_m
    bottomPelletTopPelletContactOperator->copyMasterToSlave(columnSolVec);
    bottomPelletCladContactOperator->copyMasterToSlave(columnSolVec);
    topPelletCladContactOperator->copyMasterToSlave(columnSolVec);


    linearSolver->solve(columnRhsVec, columnSolVec);

    // u^s = C u^m + d
    bottomPelletTopPelletContactOperator->copyMasterToSlave(columnSolVec);
    bottomPelletCladContactOperator->copyMasterToSlave(columnSolVec);
    topPelletCladContactOperator->copyMasterToSlave(columnSolVec);

    bottomPelletTopPelletContactOperator->addShiftToSlave(columnSolVec);
    bottomPelletCladContactOperator->addShiftToSlave(columnSolVec);
    topPelletCladContactOperator->addShiftToSlave(columnSolVec);

    computeStressTensor(meshAdapter, columnSolVec, 
        sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy,
        sigma_eff, 1.0e6, 0.3,
        referenceTemperature, thermalExpansionCoefficient, tempVec);

    activeSetVec->setToScalar(-1.0);

    std::vector<AMP::Mesh::MeshElementID> const * bottomPelletTopPelletPointerToActiveSet;
    bottomPelletTopPelletContactOperator->getActiveSet(bottomPelletTopPelletPointerToActiveSet);
    size_t const bottomPelletTopPelletSizeOfActiveSetBeforeUpdate = bottomPelletTopPelletPointerToActiveSet->size();

    std::vector<size_t> bottomPelletTopPelletActiveSetTempDOFsIndicesBeforeUpdate;
    tempDofManager->getDOFs(*bottomPelletTopPelletPointerToActiveSet, bottomPelletTopPelletActiveSetTempDOFsIndicesBeforeUpdate);
    AMP_ASSERT( bottomPelletTopPelletActiveSetTempDOFsIndicesBeforeUpdate.size() == bottomPelletTopPelletSizeOfActiveSetBeforeUpdate );
    std::vector<double> bottomPelletTopPelletValuesForActiveSet(bottomPelletTopPelletSizeOfActiveSetBeforeUpdate, 2.0); 
    activeSetVec->setLocalValuesByGlobalID(bottomPelletTopPelletSizeOfActiveSetBeforeUpdate, &(bottomPelletTopPelletActiveSetTempDOFsIndicesBeforeUpdate[0]), &(bottomPelletTopPelletValuesForActiveSet[0]));

    std::vector<size_t> bottomPelletTopPelletActiveSetDispDOFsIndicesBeforeUpdate;
    dispDofManager->getDOFs(*bottomPelletTopPelletPointerToActiveSet, bottomPelletTopPelletActiveSetDispDOFsIndicesBeforeUpdate);
    AMP_ASSERT( bottomPelletTopPelletActiveSetDispDOFsIndicesBeforeUpdate.size() == 3*bottomPelletTopPelletSizeOfActiveSetBeforeUpdate );
//
    std::vector<AMP::Mesh::MeshElementID> const * bottomPelletCladPointerToActiveSet;
    bottomPelletCladContactOperator->getActiveSet(bottomPelletCladPointerToActiveSet);
    size_t const bottomPelletCladSizeOfActiveSetBeforeUpdate = bottomPelletCladPointerToActiveSet->size();

    std::vector<size_t> bottomPelletCladActiveSetTempDOFsIndicesBeforeUpdate;
    tempDofManager->getDOFs(*bottomPelletCladPointerToActiveSet, bottomPelletCladActiveSetTempDOFsIndicesBeforeUpdate);
    AMP_ASSERT( bottomPelletCladActiveSetTempDOFsIndicesBeforeUpdate.size() == bottomPelletCladSizeOfActiveSetBeforeUpdate );
    std::vector<double> bottomPelletCladValuesForActiveSet(bottomPelletCladSizeOfActiveSetBeforeUpdate, 2.0); 
    activeSetVec->setLocalValuesByGlobalID(bottomPelletCladSizeOfActiveSetBeforeUpdate, &(bottomPelletCladActiveSetTempDOFsIndicesBeforeUpdate[0]), &(bottomPelletCladValuesForActiveSet[0]));

    std::vector<size_t> bottomPelletCladActiveSetDispDOFsIndicesBeforeUpdate;
    dispDofManager->getDOFs(*bottomPelletCladPointerToActiveSet, bottomPelletCladActiveSetDispDOFsIndicesBeforeUpdate);
    AMP_ASSERT( bottomPelletCladActiveSetDispDOFsIndicesBeforeUpdate.size() == 3*bottomPelletCladSizeOfActiveSetBeforeUpdate );
//
    std::vector<AMP::Mesh::MeshElementID> const * topPelletCladPointerToActiveSet;
    topPelletCladContactOperator->getActiveSet(topPelletCladPointerToActiveSet);
    size_t const topPelletCladSizeOfActiveSetBeforeUpdate = topPelletCladPointerToActiveSet->size();

    std::vector<size_t> topPelletCladActiveSetTempDOFsIndicesBeforeUpdate;
    tempDofManager->getDOFs(*topPelletCladPointerToActiveSet, topPelletCladActiveSetTempDOFsIndicesBeforeUpdate);
    AMP_ASSERT( topPelletCladActiveSetTempDOFsIndicesBeforeUpdate.size() == topPelletCladSizeOfActiveSetBeforeUpdate );
    std::vector<double> topPelletCladValuesForActiveSet(topPelletCladSizeOfActiveSetBeforeUpdate, 2.0); 
    activeSetVec->setLocalValuesByGlobalID(topPelletCladSizeOfActiveSetBeforeUpdate, &(topPelletCladActiveSetTempDOFsIndicesBeforeUpdate[0]), &(topPelletCladValuesForActiveSet[0]));

    std::vector<size_t> topPelletCladActiveSetDispDOFsIndicesBeforeUpdate;
    dispDofManager->getDOFs(*topPelletCladPointerToActiveSet, topPelletCladActiveSetDispDOFsIndicesBeforeUpdate);
    AMP_ASSERT( topPelletCladActiveSetDispDOFsIndicesBeforeUpdate.size() == 3*topPelletCladSizeOfActiveSetBeforeUpdate );


#ifdef USE_EXT_SILO
{
    meshAdapter->displaceMesh(columnSolVec);
    char outFileName[256];
    sprintf(outFileName, "TATA_%d", 0);
//    siloWriter->writeFile(outFileName, (activeSetIteration+1)+(thermalLoadingIteration)*maxActiveSetIterations);
    siloWriter->writeFile(outFileName, TOTO_count);
    columnSolVec->scale(-1.0);
    meshAdapter->displaceMesh(columnSolVec);
    columnSolVec->scale(-1.0);
}
#endif

    size_t nChangesInActiveSet = 0;
    nChangesInActiveSet += bottomPelletTopPelletContactOperator->updateActiveSet(columnSolVec);
    nChangesInActiveSet += bottomPelletCladContactOperator->updateActiveSet(columnSolVec);
    nChangesInActiveSet += topPelletCladContactOperator->updateActiveSet(columnSolVec);

    suckItVec->zero();
    surfaceTractionVec->zero();
    normalVectorVec->zero();

    size_t const bottomPelletTopPelletSizeOfActiveSetAfterUpdate = bottomPelletTopPelletPointerToActiveSet->size();

    std::vector<double> const * bottomPelletTopPelletSlaveVerticesNormalVector;
    std::vector<double> const * bottomPelletTopPelletSlaveVerticesSurfaceTraction;
    bottomPelletTopPelletContactOperator->getSlaveVerticesNormalVectorAndSurfaceTraction(bottomPelletTopPelletSlaveVerticesNormalVector, bottomPelletTopPelletSlaveVerticesSurfaceTraction);
    AMP_ASSERT( bottomPelletTopPelletSlaveVerticesSurfaceTraction->size() == 3*bottomPelletTopPelletSizeOfActiveSetBeforeUpdate);
    AMP_ASSERT( bottomPelletTopPelletSlaveVerticesNormalVector->size() == 3*bottomPelletTopPelletSizeOfActiveSetAfterUpdate);
    surfaceTractionVec->setLocalValuesByGlobalID(3*bottomPelletTopPelletSizeOfActiveSetBeforeUpdate, &(bottomPelletTopPelletActiveSetDispDOFsIndicesBeforeUpdate[0]), &((*bottomPelletTopPelletSlaveVerticesSurfaceTraction)[0]));
    normalVectorVec->setLocalValuesByGlobalID(3*bottomPelletTopPelletSizeOfActiveSetBeforeUpdate, &(bottomPelletTopPelletActiveSetDispDOFsIndicesBeforeUpdate[0]), &((*bottomPelletTopPelletSlaveVerticesNormalVector)[0]));

    std::vector<double> bottomPelletTopPelletSurfaceTractionDOTnormalVector(bottomPelletTopPelletSizeOfActiveSetBeforeUpdate);
    for (size_t kk = 0; kk < bottomPelletTopPelletSizeOfActiveSetBeforeUpdate; ++kk) {
      bottomPelletTopPelletSurfaceTractionDOTnormalVector[kk] = - compute_scalar_product(&((*bottomPelletTopPelletSlaveVerticesSurfaceTraction)[3*kk]), &((*bottomPelletTopPelletSlaveVerticesNormalVector)[3*kk]));
    } // end for kk
    suckItVec->setLocalValuesByGlobalID(bottomPelletTopPelletSizeOfActiveSetBeforeUpdate, &(bottomPelletTopPelletActiveSetTempDOFsIndicesBeforeUpdate[0]), &(bottomPelletTopPelletSurfaceTractionDOTnormalVector[0]));
//
    size_t const bottomPelletCladSizeOfActiveSetAfterUpdate = bottomPelletCladPointerToActiveSet->size();

    std::vector<double> const * bottomPelletCladSlaveVerticesNormalVector;
    std::vector<double> const * bottomPelletCladSlaveVerticesSurfaceTraction;
    bottomPelletCladContactOperator->getSlaveVerticesNormalVectorAndSurfaceTraction(bottomPelletCladSlaveVerticesNormalVector, bottomPelletCladSlaveVerticesSurfaceTraction);
    AMP_ASSERT( bottomPelletCladSlaveVerticesSurfaceTraction->size() == 3*bottomPelletCladSizeOfActiveSetBeforeUpdate);
    AMP_ASSERT( bottomPelletCladSlaveVerticesNormalVector->size() == 3*bottomPelletCladSizeOfActiveSetAfterUpdate);
    surfaceTractionVec->setLocalValuesByGlobalID(3*bottomPelletCladSizeOfActiveSetBeforeUpdate, &(bottomPelletCladActiveSetDispDOFsIndicesBeforeUpdate[0]), &((*bottomPelletCladSlaveVerticesSurfaceTraction)[0]));
    normalVectorVec->setLocalValuesByGlobalID(3*bottomPelletCladSizeOfActiveSetBeforeUpdate, &(bottomPelletCladActiveSetDispDOFsIndicesBeforeUpdate[0]), &((*bottomPelletCladSlaveVerticesNormalVector)[0]));

    std::vector<double> bottomPelletCladSurfaceTractionDOTnormalVector(bottomPelletCladSizeOfActiveSetBeforeUpdate);
    for (size_t kk = 0; kk < bottomPelletCladSizeOfActiveSetBeforeUpdate; ++kk) {
      bottomPelletCladSurfaceTractionDOTnormalVector[kk] = - compute_scalar_product(&((*bottomPelletCladSlaveVerticesSurfaceTraction)[3*kk]), &((*bottomPelletCladSlaveVerticesNormalVector)[3*kk]));
    } // end for kk
    suckItVec->setLocalValuesByGlobalID(bottomPelletCladSizeOfActiveSetBeforeUpdate, &(bottomPelletCladActiveSetTempDOFsIndicesBeforeUpdate[0]), &(bottomPelletCladSurfaceTractionDOTnormalVector[0]));
//
    size_t const topPelletCladSizeOfActiveSetAfterUpdate = topPelletCladPointerToActiveSet->size();

    std::vector<double> const * topPelletCladSlaveVerticesNormalVector;
    std::vector<double> const * topPelletCladSlaveVerticesSurfaceTraction;
    topPelletCladContactOperator->getSlaveVerticesNormalVectorAndSurfaceTraction(topPelletCladSlaveVerticesNormalVector, topPelletCladSlaveVerticesSurfaceTraction);
    AMP_ASSERT( topPelletCladSlaveVerticesSurfaceTraction->size() == 3*topPelletCladSizeOfActiveSetBeforeUpdate);
    AMP_ASSERT( topPelletCladSlaveVerticesNormalVector->size() == 3*topPelletCladSizeOfActiveSetAfterUpdate);
    surfaceTractionVec->setLocalValuesByGlobalID(3*topPelletCladSizeOfActiveSetBeforeUpdate, &(topPelletCladActiveSetDispDOFsIndicesBeforeUpdate[0]), &((*topPelletCladSlaveVerticesSurfaceTraction)[0]));
    normalVectorVec->setLocalValuesByGlobalID(3*topPelletCladSizeOfActiveSetBeforeUpdate, &(topPelletCladActiveSetDispDOFsIndicesBeforeUpdate[0]), &((*topPelletCladSlaveVerticesNormalVector)[0]));

    std::vector<double> topPelletCladSurfaceTractionDOTnormalVector(topPelletCladSizeOfActiveSetBeforeUpdate);
    for (size_t kk = 0; kk < topPelletCladSizeOfActiveSetBeforeUpdate; ++kk) {
      topPelletCladSurfaceTractionDOTnormalVector[kk] = - compute_scalar_product(&((*topPelletCladSlaveVerticesSurfaceTraction)[3*kk]), &((*topPelletCladSlaveVerticesNormalVector)[3*kk]));
    } // end for kk
    suckItVec->setLocalValuesByGlobalID(topPelletCladSizeOfActiveSetBeforeUpdate, &(topPelletCladActiveSetTempDOFsIndicesBeforeUpdate[0]), &(topPelletCladSurfaceTractionDOTnormalVector[0]));
    

oldSolVec->subtract(columnSolVec, oldSolVec);
#ifdef USE_EXT_SILO
    meshAdapter->displaceMesh(columnSolVec);
    char outFileName[256];
    sprintf(outFileName, "TATA_%d", 0);
//    siloWriter->writeFile(outFileName, (activeSetIteration+1)+(thermalLoadingIteration)*maxActiveSetIterations);
    siloWriter->writeFile(outFileName, TOTO_count);
    columnSolVec->scale(-1.0);
    meshAdapter->displaceMesh(columnSolVec);
    columnSolVec->scale(-1.0);
#endif

    if (!rank) { std::cout<<nChangesInActiveSet<<" CHANGES IN ACTIVE SET\n"; }

double errL1Norm = oldSolVec->L1Norm();
double solL1Norm = columnSolVec->L1Norm();
double relErrL1Norm = errL1Norm / solL1Norm;
double errL2Norm = oldSolVec->L2Norm();
double solL2Norm = columnSolVec->L2Norm();
double relErrL2Norm = errL2Norm / solL2Norm;
    if (!rank) { std::cout<<"ERROR L1 NORM "<<errL1Norm<<" ("<<100.0*relErrL1Norm<<"%)    "; }
    if (!rank) { std::cout<<"ERROR L2 NORM "<<errL2Norm<<" ("<<100.0*relErrL2Norm<<"%)  \n"; }
oldSolVec->copyVector(columnSolVec);

    if (nChangesInActiveSet == 0) { break; }
//    AMP_ASSERT( activeSetIteration != maxActiveSetIterations - 1 );
    if ( activeSetIteration == maxActiveSetIterations - 1 ) {
      if (!rank) { std::cout<<"!!!!!! ACTIVE SET ITERATIONS DID NOT CONVERGE !!!!!!!!\n"; }
    } // end if
  } // end for

} // end for

  meshAdapter->displaceMesh(columnSolVec);

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
  exeNames.push_back("testNodeToFaceContactOperator-5");

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



