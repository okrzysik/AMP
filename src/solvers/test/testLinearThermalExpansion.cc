
#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "boost/shared_ptr.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#include "ampmesh/SiloIO.h"

#include "operators/ColumnOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/LinearBVPOperator.h"
#include "operators/mechanics/ConstructLinearMechanicsRHSVector.h"

#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/TrilinosMLSolver.h"

/** Post-processing: Move the mesh using the given displacement field */
void deformMesh(AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter,
    AMP::LinearAlgebra::Vector::shared_ptr mechSolVec) {
  AMP::Mesh::DOFMap::shared_ptr dof_map = meshAdapter->getDOFMap( mechSolVec->getVariable() );

  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator nd = meshAdapter->beginOwnedNode();
  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator end_nd = meshAdapter->endOwnedNode();

  std::vector <unsigned int> dofIds(3);
  dofIds[0] = 0; dofIds[1] = 1; dofIds[2] = 2;

  for( ; nd != end_nd; ++nd) {
    std::vector<unsigned int> ndGlobalIds;
    dof_map->getDOFs(*nd, ndGlobalIds, dofIds);

    double xDisp = mechSolVec->getValueByGlobalID(ndGlobalIds[0]);
    double yDisp = mechSolVec->getValueByGlobalID(ndGlobalIds[1]);
    double zDisp = mechSolVec->getValueByGlobalID(ndGlobalIds[2]);

    nd->translate(xDisp, yDisp, zDisp);
  }//end for nd
}

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "pellet" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator = boost::dynamic_pointer_cast<
    AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "MechanicsBVPOperator", input_db, elementPhysicsModel));

  boost::shared_ptr<AMP::Database> temperatureRhsDatabase = input_db->getDatabase("TemperatureRHS");

  AMP::LinearAlgebra::Variable::shared_ptr dispVar = bvpOperator->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr tempVar(new 
      AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>("temp", meshAdapter) );

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr mechSolVec = meshAdapter->createVector( dispVar );
  AMP::LinearAlgebra::Vector::shared_ptr mechRhsVec = meshAdapter->createVector( dispVar );
  AMP::LinearAlgebra::Vector::shared_ptr mechResVec = meshAdapter->createVector( dispVar );

  AMP::LinearAlgebra::Vector::shared_ptr currTempVec = meshAdapter->createVector( tempVar );
  AMP::LinearAlgebra::Vector::shared_ptr prevTempVec = meshAdapter->createVector( tempVar );

  mechSolVec->setToScalar(0.0);
  mechResVec->setToScalar(0.0);

  currTempVec->setToScalar(500.0);
  prevTempVec->setToScalar(300.0);

  computeTemperatureRhsVector(meshAdapter, temperatureRhsDatabase, tempVar, dispVar, currTempVec, prevTempVec, mechRhsVec);
  bvpOperator->modifyRHSvector(mechRhsVec);

  boost::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase("LinearSolver"); 
  boost::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase("Preconditioner"); 
  boost::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(new AMP::Solver::TrilinosMLSolverParameters(pcSolver_db));
  pcSolverParams->d_pOperator = bvpOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(new AMP::Solver::TrilinosMLSolver(pcSolverParams));

  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(linearSolver_db));
  linearSolverParams->d_pOperator = bvpOperator;
  linearSolverParams->d_comm = AMP_COMM_WORLD;
  linearSolverParams->d_pPreconditioner = pcSolver;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new AMP::Solver::PetscKrylovSolver(linearSolverParams));

  linearSolver->solve(mechRhsVec, mechSolVec);

#ifdef USE_SILO
  manager->registerVectorAsData ( mechSolVec, "Displacement" );
  manager->writeFile<AMP::Mesh::SiloIO> ( exeName , 1 );
  deformMesh(meshAdapter, mechSolVec);
  manager->writeFile<AMP::Mesh::SiloIO> ( exeName , 2 );
#endif

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testLinearThermalExpansion";

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


