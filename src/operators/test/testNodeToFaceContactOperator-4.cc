
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
#include "operators/mechanics/ConstructLinearMechanicsRHSVector.h"
#include "operators/contact/NodeToFaceContactOperator.h"
#include "operators/mechanics/IsotropicElasticModel.h"

#include "solvers/ColumnSolver.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/ConstraintsEliminationSolver.h"

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
  boost::shared_ptr<AMP::Operator::MechanicsModelParameters> masterMechanicsMaterialModelParams(new AMP::Operator::MechanicsModelParameters(model_db));
  boost::shared_ptr<AMP::Operator::MechanicsMaterialModel> masterMechanicsMaterialModel(new AMP::Operator::IsotropicElasticModel(masterMechanicsMaterialModelParams));

  // Build the contact operator
  AMP_INSIST(input_db->keyExists("ContactOperator"), "Key ''ContactOperator'' is missing!");
  boost::shared_ptr<AMP::Database> contact_db = input_db->getDatabase("ContactOperator");
  boost::shared_ptr<AMP::Operator::ContactOperatorParameters> contactOperatorParams(new AMP::Operator::ContactOperatorParameters(contact_db));
  contactOperatorParams->d_DOFsPerNode = dofsPerNode;
  contactOperatorParams->d_DOFManager = dispDofManager;
  contactOperatorParams->d_GlobalComm = globalComm;
  contactOperatorParams->d_Mesh = meshAdapter;
  contactOperatorParams->d_MasterMechanicsMaterialModel = masterMechanicsMaterialModel;
  contactOperatorParams->reset(); // got segfault at constructor since d_Mesh was pointing to NULL
  boost::shared_ptr<AMP::Operator::NodeToFaceContactOperator> contactOperator(new AMP::Operator::NodeToFaceContactOperator(contactOperatorParams));
  contactOperator->initialize();
  
  // Build the master and slave operators
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> masterBVPOperator;
  AMP::Mesh::MeshID masterMeshID = contactOperator->getMasterMeshID();
  AMP::Mesh::Mesh::shared_ptr masterMeshAdapter = meshAdapter->Subset(masterMeshID);
  // NB: need to rotate the mesh before building mechanics op 
//  rotateMesh(masterMeshAdapter);

  if (masterMeshAdapter.get() != NULL) {
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
    masterBVPOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(masterMeshAdapter,
                                                                                                                                     "MasterBVPOperator",
                                                                                                                                     input_db,
                                                                                                                                     masterElementPhysicsModel));
    columnOperator->append(masterBVPOperator);

    boost::shared_ptr<AMP::Database> masterSolver_db = columnPreconditioner_db->getDatabase("MasterSolver"); 
    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> masterSolverParams(new AMP::Solver::PetscKrylovSolverParameters(masterSolver_db));
    masterSolverParams->d_pOperator = masterBVPOperator;
    masterSolverParams->d_comm = masterMeshAdapter->getComm();
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> masterSolver(new AMP::Solver::PetscKrylovSolver(masterSolverParams));
    columnPreconditioner->append(masterSolver);
  } // end if

  boost::shared_ptr<AMP::Operator::LinearBVPOperator> slaveBVPOperator;

  AMP::Mesh::MeshID slaveMeshID = contactOperator->getSlaveMeshID();
  AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter = meshAdapter->Subset(slaveMeshID);
  if (slaveMeshAdapter.get() != NULL) {
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;
    slaveBVPOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter,
                                                                                                                                    "SlaveBVPOperator",
                                                                                                                                    input_db,
                                                                                                                                    slaveElementPhysicsModel));
    columnOperator->append(slaveBVPOperator);

    boost::shared_ptr<AMP::Database> slaveSolver_db = columnPreconditioner_db->getDatabase("SlaveSolver"); 
    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> slaveSolverParams(new AMP::Solver::PetscKrylovSolverParameters(slaveSolver_db));
    slaveSolverParams->d_pOperator = slaveBVPOperator;
    slaveSolverParams->d_comm = slaveMeshAdapter->getComm();
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> slaveSolver(new AMP::Solver::PetscKrylovSolver(slaveSolverParams));
    columnPreconditioner->append(slaveSolver);
  } // end if

  boost::shared_ptr<AMP::Database> contactPreconditioner_db = columnPreconditioner_db->getDatabase("ContactPreconditioner"); 
  boost::shared_ptr<AMP::Solver::ConstraintsEliminationSolverParameters> contactPreconditionerParams(new 
      AMP::Solver::ConstraintsEliminationSolverParameters(contactPreconditioner_db));
  contactPreconditionerParams->d_pOperator = contactOperator;
  boost::shared_ptr<AMP::Solver::ConstraintsEliminationSolver> contactPreconditioner(new AMP::Solver::ConstraintsEliminationSolver(contactPreconditionerParams));
  columnPreconditioner->append(contactPreconditioner);

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

  tempVec->setToScalar(500.0);
//AMP::Mesh::MeshIterator meshIterator = meshAdapter->getIterator(AMP::Mesh::Vertex),
std::vector<double> vertexCoord;
std::vector<size_t> DOFsIndices;
double temperatureFuelOuterRadius = input_db->getDouble("TemperatureFuelOuterRadius"); 
double temperatureCladInnerRadius = input_db->getDouble("TemperatureCladInnerRadius"); 
double temperatureCladOuterRadius = input_db->getDouble("TemperatureCladOuterRadius"); 
double heatGenerationRate = input_db->getDouble("HeatGenerationRate");
double fuelOuterRadius = input_db->getDouble("FuelOuterRadius");
double cladInnerRadius = input_db->getDouble("CladInnerRadius");
double cladOuterRadius = input_db->getDouble("CladOuterRadius");
double fuelOuterRadiusSquared = fuelOuterRadius * fuelOuterRadius;
double fuelThermalConductivity = input_db->getDouble("FuelThermalConductivity");
double temperatureCenterLine = temperatureFuelOuterRadius + heatGenerationRate * fuelOuterRadiusSquared / (4.0 * fuelThermalConductivity);
double referenceTemperature = input_db->getDouble("ReferenceTemperature");
if (!rank) { std::cout<<"temperatureCenterLine="<<temperatureCenterLine<<"\n"; }
  refTempVec->setToScalar(referenceTemperature);
  tempVec->setToScalar(referenceTemperature);

{ // temperature in fuel
AMP::Mesh::MeshIterator meshIterator = slaveMeshAdapter->getIterator(AMP::Mesh::Vertex);
AMP::Mesh::MeshIterator meshIterator_begin = meshIterator.begin();
AMP::Mesh::MeshIterator meshIterator_end = meshIterator.end();
for (meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator) {
  vertexCoord = meshIterator->coord();
  double radiusSquared = vertexCoord[0]*vertexCoord[0] + vertexCoord[1]*vertexCoord[1];
  double temperature = temperatureCenterLine - heatGenerationRate * radiusSquared / (4.0 * fuelThermalConductivity);
  tempDofManager->getDOFs(meshIterator->globalID(), DOFsIndices);
  AMP_ASSERT(DOFsIndices.size() == 1);
  tempVec->setLocalValuesByGlobalID(1, &(DOFsIndices[0]), &temperature);
} // end for
}

{ // temperature in clad
AMP::Mesh::MeshIterator meshIterator = masterMeshAdapter->getIterator(AMP::Mesh::Vertex);
AMP::Mesh::MeshIterator meshIterator_begin = meshIterator.begin();
AMP::Mesh::MeshIterator meshIterator_end = meshIterator.end();
for (meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator) {
  vertexCoord = meshIterator->coord();
  double radius = std::sqrt(vertexCoord[0]*vertexCoord[0] + vertexCoord[1]*vertexCoord[1]);
  AMP_ASSERT((radius >= cladInnerRadius) && (radius <= cladOuterRadius));
  double temperature = temperatureCladInnerRadius + (temperatureCladOuterRadius - temperatureCladInnerRadius) * (radius - cladInnerRadius) / (cladOuterRadius - cladOuterRadius);
  tempDofManager->getDOFs(meshIterator->globalID(), DOFsIndices);
  AMP_ASSERT(DOFsIndices.size() == 1);
  tempVec->setLocalValuesByGlobalID(1, &(DOFsIndices[0]), &temperature);
} // end for
}

AMP::LinearAlgebra::VS_Mesh slaveVectorSelector(slaveMeshAdapter);
AMP::LinearAlgebra::Vector::shared_ptr slaveTempVec = tempVec->select(slaveVectorSelector, tempVar->getName());
//slaveTempVec->setToScalar(900.0);
boost::shared_ptr<AMP::Database> tmp_db = temperatureRhs_db->getDatabase("RhsMaterialModel");
double thermalExpansionCoefficient = tmp_db->getDouble("THERMAL_EXPANSION_COEFFICIENT");
  contactOperator->uglyHack(tempVec, tempDofManager, thermalExpansionCoefficient, referenceTemperature);

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
  contactOperator->updateActiveSet(nullVec, skipDisplaceMesh);

  
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
    sprintf(outFileName, "TOTO_%d", 0);
    siloWriter->writeFile(outFileName, 0);
  }
#endif
  oldSolVec->copyVector(columnSolVec);

  columnSolVec->zero();
  columnOperator->append(contactOperator);

  // Build a matrix shell operator to use the column operator with the petsc krylov solvers
  boost::shared_ptr<AMP::Database> matrixShellDatabase = input_db->getDatabase("MatrixShellOperator");
  boost::shared_ptr<AMP::Operator::OperatorParameters> matrixShellParams(new AMP::Operator::OperatorParameters(matrixShellDatabase));
  boost::shared_ptr<AMP::Operator::PetscMatrixShellOperator> matrixShellOperator(new AMP::Operator::PetscMatrixShellOperator(matrixShellParams));

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

  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = matrixShellOperator;
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pPreconditioner = columnPreconditioner;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));
//  linearSolver->setZeroInitialGuess(true);
  linearSolver->setInitialGuess(columnSolVec);

  AMP::LinearAlgebra::Vector::shared_ptr fullThermalLoadingTempMinusRefTempVec = tempVec->cloneVector();
  fullThermalLoadingTempMinusRefTempVec->subtract(tempVec, refTempVec);

size_t const maxThermalLoadingIterations = input_db->getIntegerWithDefault("maxThermalLoadingIterations", 5);
for (size_t thermalLoadingIteration = 0; thermalLoadingIteration < maxThermalLoadingIterations; ++thermalLoadingIteration) {
  if (!rank) { std::cout<<"THERMAL LOADING "<<thermalLoadingIteration+1<<"/"<<maxThermalLoadingIterations<<"\n"; }
  double scalingFactor = static_cast<double>(thermalLoadingIteration+1) / static_cast<double>(maxThermalLoadingIterations);
  tempVec->axpy(scalingFactor, fullThermalLoadingTempMinusRefTempVec, refTempVec);

  size_t const maxActiveSetIterations = input_db->getIntegerWithDefault("maxActiveSetIterations", 5);
  for (size_t activeSetIteration = 0; activeSetIteration < maxActiveSetIterations; ++activeSetIteration) {
    if (!rank) { std::cout<<"ACTIVE SET ITERATION #"<<activeSetIteration+1<<std::endl; }

    columnSolVec->zero();
    columnRhsVec->zero();

    // compute thermal load f
    computeTemperatureRhsVector(meshAdapter, temperatureRhs_db, tempVar, dispVar, tempVec, refTempVec, columnRhsVec);

    // apply dirichlet rhs correction on f
    if (masterBVPOperator.get() != NULL) {
      masterBVPOperator->modifyRHSvector(columnRhsVec);
    } // end if
    if (slaveBVPOperator.get() != NULL) {
      slaveBVPOperator->modifyRHSvector(columnRhsVec);
    } // end if

    // get d
//    AMP::LinearAlgebra::Vector::shared_ptr contactShiftVec = createVector(dispDofManager, columnVar, split);
    contactShiftVec->zero();
    contactOperator->addShiftToSlave(contactShiftVec);
//  contactOperator->addShiftToSlave(columnSolVec);

    // compute - Kd
    AMP::LinearAlgebra::Vector::shared_ptr rhsCorrectionVec = createVector(dispDofManager, columnVar, split);
    rhsCorrectionVec->zero();
    masterBVPOperator->apply(nullVec, contactShiftVec, rhsCorrectionVec, -1.0, 0.0);
    slaveBVPOperator->apply(nullVec, contactShiftVec, rhsCorrectionVec, -1.0, 0.0);
//  columnOperator->apply(nullVec, columnSolVec, rhsCorrectionVec, -1.0, 0.0);
//  columnOperator->append(contactOperator);

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

    computeStressTensor(meshAdapter, columnSolVec, 
        sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy,
        sigma_eff, 1.0e6, 0.3,
        referenceTemperature, thermalExpansionCoefficient, tempVec);

    std::vector<AMP::Mesh::MeshElementID> const * pointerToActiveSet;
    contactOperator->getActiveSet(pointerToActiveSet);
    size_t const sizeOfActiveSetBeforeUpdate = pointerToActiveSet->size();

    std::vector<size_t> activeSetTempDOFsIndicesBeforeUpdate;
    tempDofManager->getDOFs(*pointerToActiveSet, activeSetTempDOFsIndicesBeforeUpdate);
    AMP_ASSERT( activeSetTempDOFsIndicesBeforeUpdate.size() == sizeOfActiveSetBeforeUpdate );
    std::vector<double> valuesForActiveSet(sizeOfActiveSetBeforeUpdate, 2.0); 
    activeSetVec->setToScalar(-1.0);
    activeSetVec->setLocalValuesByGlobalID(sizeOfActiveSetBeforeUpdate, &(activeSetTempDOFsIndicesBeforeUpdate[0]), &(valuesForActiveSet[0]));

    std::vector<size_t> activeSetDispDOFsIndicesBeforeUpdate;
    dispDofManager->getDOFs(*pointerToActiveSet, activeSetDispDOFsIndicesBeforeUpdate);
    AMP_ASSERT( activeSetDispDOFsIndicesBeforeUpdate.size() == 3*sizeOfActiveSetBeforeUpdate );
    
#ifdef USE_EXT_SILO
{
    meshAdapter->displaceMesh(columnSolVec);
    char outFileName[256];
    sprintf(outFileName, "TOTO_%d", 0);
    siloWriter->writeFile(outFileName, (activeSetIteration+1)+(thermalLoadingIteration)*maxActiveSetIterations);
    columnSolVec->scale(-1.0);
    meshAdapter->displaceMesh(columnSolVec);
    columnSolVec->scale(-1.0);
}
#endif

    size_t nChangesInActiveSet = contactOperator->updateActiveSet(columnSolVec);

    size_t const sizeOfActiveSetAfterUpdate = pointerToActiveSet->size();
//    std::vector<size_t> activeSetDOFsIndicesAfterUpdate;
//    tempDofManager->getDOFs(*pointerToActiveSet, activeSetDOFsIndicesAfterUpdate);
//    AMP_ASSERT( activeSetDOFsIndicesAfterUpdate.size() == sizeOfActiveSetAfterUpdate );
//    std::vector<double> valuesForActiveSet(pointerToActiveSet->size(), 2.0); 
//    activeSetVec->setToScalar(-1.0);
//    activeSetVec->setLocalValuesByGlobalID(sizeOfActiveSetAfterUpdate, &(activeSetDOFsIndicesAfterUpdate[0]), &(valuesForActiveSet[0]));

    std::vector<double> const * slaveVerticesNormalVector;
    std::vector<double> const * slaveVerticesSurfaceTraction;
    contactOperator->getSlaveVerticesNormalVectorAndSurfaceTraction(slaveVerticesNormalVector, slaveVerticesSurfaceTraction);
    AMP_ASSERT( slaveVerticesSurfaceTraction->size() == 3*sizeOfActiveSetBeforeUpdate);
    AMP_ASSERT( slaveVerticesNormalVector->size() == 3*sizeOfActiveSetAfterUpdate);
    surfaceTractionVec->zero();
    surfaceTractionVec->setLocalValuesByGlobalID(3*sizeOfActiveSetBeforeUpdate, &(activeSetDispDOFsIndicesBeforeUpdate[0]), &((*slaveVerticesSurfaceTraction)[0]));
    normalVectorVec->zero();
    normalVectorVec->setLocalValuesByGlobalID(3*sizeOfActiveSetBeforeUpdate, &(activeSetDispDOFsIndicesBeforeUpdate[0]), &((*slaveVerticesNormalVector)[0]));

    std::vector<double> surfaceTractionDOTnormalVector(sizeOfActiveSetBeforeUpdate);
    for (size_t kk = 0; kk < sizeOfActiveSetBeforeUpdate; ++kk) {
      surfaceTractionDOTnormalVector[kk] = - compute_scalar_product(&((*slaveVerticesSurfaceTraction)[3*kk]), &((*slaveVerticesNormalVector)[3*kk]));
    } // end for kk
    suckItVec->zero();
    suckItVec->setLocalValuesByGlobalID(sizeOfActiveSetBeforeUpdate, &(activeSetTempDOFsIndicesBeforeUpdate[0]), &(surfaceTractionDOTnormalVector[0]));
    
    
//why_cant_we_be_friend(masterMeshAdapter, columnSolVec);

oldSolVec->subtract(columnSolVec, oldSolVec);
#ifdef USE_EXT_SILO
    meshAdapter->displaceMesh(columnSolVec);
/*    siloWriter->registerVector(columnSolVec, meshAdapter, AMP::Mesh::Vertex, "Solution");
    siloWriter->registerVector(tempVec, meshAdapter, AMP::Mesh::Vertex, "Temperature");
    siloWriter->registerVector(sigma_eff, meshAdapter, AMP::Mesh::Vertex, "vonMises");
    siloWriter->registerVector(activeSetVec, meshAdapter, AMP::Mesh::Vertex, "Contact");
    siloWriter->registerVector(oldSolVec, meshAdapter, AMP::Mesh::Vertex, "Error");
    siloWriter->registerVector(surfaceTractionVec, meshAdapter, AMP::Mesh::Vertex, "Traction");
*/
    char outFileName[256];
    sprintf(outFileName, "TOTO_%d", 0);
    siloWriter->writeFile(outFileName, (activeSetIteration+1)+(thermalLoadingIteration)*maxActiveSetIterations);
    columnSolVec->scale(-1.0);
    meshAdapter->displaceMesh(columnSolVec);
    columnSolVec->scale(-1.0);
#endif
//    for (std::vector<AMP::Mesh::MeshElementID>::iterator activeSetIterator = pointerToActiveSet->begin(); activeSetIterator != pointerToActiveSet->end(); ++activeSetIterator) {
//    } // end for
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
  exeNames.push_back("testNodeToFaceContactOperator-4");

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



