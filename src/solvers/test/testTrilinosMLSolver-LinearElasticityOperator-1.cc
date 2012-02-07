
#include <iostream>
#include <string>

#include <cassert>

#include <fstream>

#include <sys/stat.h>

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

/* Boost files */
#include "boost/shared_ptr.hpp"

/* libMesh files */
#include "mesh.h"
#include "mesh_generation.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "petsc_matrix.h"
#include "petsc_vector.h"
#include "dense_matrix.h"
#include "linear_implicit_system.h"
#include "elem.h"

/* AMP files */
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "ampmesh/MeshVariable.h"
#include "materials/Material.h"


#include "operators/mechanics/IsotropicElasticModel.h"
#include "operators/mechanics/MechanicsLinearElement.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "vectors/Vector.h"
#include "ampmesh/SiloIO.h"


#include "../TrilinosMLSolver.h"


void linearElasticTest(AMP::UnitTest *ut )
{
  std::string exeName("testTrilinosMLSolver-LinearElasticityOperator-1");
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( mgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "bar" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														 "MechanicsBVPOperator",
														 input_db,
														 elementPhysicsModel));

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
															 "Load_Boundary",
															 input_db,
															 dummyModel));
  //This has an in-place apply. So, it has an empty input variable and
  //the output variable is the same as what it is operating on. 
  dirichletVecOp->setVariable(bvpOperator->getOutputVariable());

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  AMP::LinearAlgebra::Vector::shared_ptr mechSolVec = meshAdapter->createVector( bvpOperator->getInputVariable() );
  AMP::LinearAlgebra::Vector::shared_ptr mechRhsVec = meshAdapter->createVector( bvpOperator->getOutputVariable() );
  AMP::LinearAlgebra::Vector::shared_ptr mechResVec = meshAdapter->createVector( bvpOperator->getOutputVariable() );

  mechSolVec->setToScalar(0.5);
  mechRhsVec->setToScalar(0.0);
  mechResVec->setToScalar(0.0);

  dirichletVecOp->apply(nullVec, nullVec, mechRhsVec, 1.0, 0.0);

  double rhsNorm = mechRhsVec->L2Norm();

  std::cout<<"RHS Norm: "<<rhsNorm<<std::endl;

  double initSolNorm = mechSolVec->L2Norm();

  std::cout<<"Initial Solution Norm: "<<initSolNorm<<std::endl;

  bvpOperator->apply(mechRhsVec, mechSolVec, mechResVec, 1.0, -1.0);

  double initResidualNorm = mechResVec->L2Norm();

  std::cout<<"Initial Residual Norm: "<<initResidualNorm<<std::endl;

  boost::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase("LinearSolver"); 

  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams(new AMP::Solver::SolverStrategyParameters(mlSolver_db));

  mlSolverParams->d_pOperator = bvpOperator;

  // create the ML solver interface
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(new AMP::Solver::TrilinosMLSolver(mlSolverParams));

  mlSolver->setZeroInitialGuess(false);

  mlSolver->solve(mechRhsVec, mechSolVec);

  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
  if( globalComm.getSize() == 1 ) {
#ifdef USE_SILO
    meshAdapter->registerVectorAsData ( mechSolVec);
    manager->writeFile<AMP::Mesh::SiloIO> ( exeName, 0 );
#endif
  }

  bvpOperator->apply(mechRhsVec, mechSolVec, mechResVec);

  double finalResidualNorm = mechResVec->L2Norm();

  std::cout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

  if(finalResidualNorm > (1e-10*initResidualNorm)) {
    ut->failure("TrilinosMLSolver successfully solves a linear elasticity problem");
  } else {
    ut->passes("TrilinosMLSolver successfully solves a linear elasticity problem");
  }

  input_db.reset();

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
        linearElasticTest(&ut);
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



