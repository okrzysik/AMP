#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "operators/VolumeIntegralOperator.h"
#include "operators/NeutronicsRhs.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "materials/Material.h"


#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshVariable.h"
#include "ampmesh/SiloIO.h"


#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"

#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/diffusion/FickSoretNonlinearFEOperator.h"
#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/ColumnOperator.h"
#include "operators/OperatorBuilder.h"

#include "../ColumnSolver.h"
#include "../PetscKrylovSolverParameters.h"
#include "../PetscKrylovSolver.h"
#include "../PetscSNESSolverParameters.h"
#include "../PetscSNESSolver.h"

#include "../TrilinosMLSolver.h"


void fickSoretTest(AMP::UnitTest *ut, std::string exeName, std::vector<double> &results)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("NumberOfMeshes"), "Key does not exist");
  int numMeshes = input_db->getInteger("NumberOfMeshes");
  
  AMP::pout<<"Num meshes = "<<numMeshes<<std::endl;

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "cylinder" );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear Fick-Soret diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearFickSoretBVPOperator"), "key missing!" );

  // Create nonlinear FickSoret BVP operator and access volume nonlinear FickSoret operator
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::Operator> nlinBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
						   "testNonlinearFickSoretBVPOperator",
						   input_db,
						   elementPhysicsModel);
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nlinBVPOp =
		  boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(nlinBVPOperator);
  boost::shared_ptr<AMP::Operator::FickSoretNonlinearFEOperator> nlinOp =
		 boost::dynamic_pointer_cast<AMP::Operator::FickSoretNonlinearFEOperator>(nlinBVPOp->getVolumeOperator());
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickOp =
		 boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinOp->getFickOperator());
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> soretOp =
		 boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinOp->getSoretOperator());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // use the linear BVP operator to create a Fick linear operator with bc's
  AMP_INSIST( input_db->keyExists("testLinearFickBVPOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::Operator> linBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
						   "testLinearFickBVPOperator",
						   input_db,
						   elementPhysicsModel);
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linBVPOp =
		  boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(linBVPOperator);
  //boost::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> linOp =
 //		 boost::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>(linBVPOp->getVolumeOperator());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // Set up input and output variables
  AMP::LinearAlgebra::Variable::shared_ptr tVar(fickOp->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
  AMP::LinearAlgebra::Variable::shared_ptr cVar(fickOp->getInputVariable(AMP::Operator::Diffusion::CONCENTRATION));
  boost::shared_ptr<AMP::LinearAlgebra::Variable> fsOutVar(nlinBVPOp->getOutputVariable());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create solution, rhs, and residual vectors

  AMP::LinearAlgebra::Vector::shared_ptr solVec = meshAdapter->createVector( cVar );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = meshAdapter->createVector( fsOutVar );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = meshAdapter->createVector( fsOutVar );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create parameters for reset test and reset fick and soret operators

  AMP::LinearAlgebra::Vector::shared_ptr tVec = meshAdapter->createVector( tVar );
  
  fickOp->setVector(0, tVec);
  soretOp->setVector(0, tVec);

  std::vector<AMP::LinearAlgebra::Vector::shared_ptr> fickFrozen  = fickOp-> getFrozen();
  std::vector<AMP::LinearAlgebra::Vector::shared_ptr> soretFrozen = soretOp->getFrozen();

  double lenscale = input_db->getDouble("LengthScale");
  soretFrozen[AMP::Operator::Diffusion::TEMPERATURE]->setToScalar(300.);  // Fill in manufactured solution
  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator iterator = meshAdapter->beginOwnedNode();
  for( ; iterator != meshAdapter->endOwnedNode(); iterator++ ) {
	double x, y, z;
	std::valarray<double> poly(10);
	x = iterator->x();
	y = iterator->y();
	z = iterator->z();
	size_t gid = iterator->globalID();
	double value = 300. + 450*(1.-(x*x/lenscale/lenscale+y*y/lenscale/lenscale));
	fickFrozen[AMP::Operator::Diffusion::TEMPERATURE]->setValueByGlobalID(gid, value);
	soretFrozen[AMP::Operator::Diffusion::TEMPERATURE]->setValueByGlobalID(gid, value);
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register some variables for plotting
  manager->registerVectorAsData ( solVec, "Solution" );
  manager->registerVectorAsData ( resVec, "Residual" );
  manager->registerVectorAsData ( fickFrozen[AMP::Operator::Diffusion::TEMPERATURE], "Temperature" );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //Initial guess

  double initialValue = input_db->getDouble("InitialValue");
  solVec->setToScalar(initialValue);
  double initialGuessNorm  = solVec->L2Norm();
  std::cout << "initial guess norm = " << initialGuessNorm <<"\n";

  nlinBVPOp->modifyInitialSolutionVector(solVec);

  initialGuessNorm  = solVec->L2Norm();
  std::cout << "initial guess norm  after apply = " << initialGuessNorm <<"\n";

  rhsVec->setToScalar(0.0);
  nlinBVPOp->modifyRHSvector(rhsVec);

  //----------------------------------------------------------------------------------------------------------------------------------------------/

  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the nonlinear solver
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nlinBVPOp;
  nonlinearSolverParams->d_pInitialGuess = solVec;

  boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  boost::shared_ptr<AMP::Database> fickPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> fickPreconditionerParams(new AMP::Solver::SolverStrategyParameters(fickPreconditioner_db));
  fickPreconditionerParams->d_pOperator = linBVPOp;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFickPreconditioner(new AMP::Solver::TrilinosMLSolver(fickPreconditionerParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register the preconditioner with the Jacobian free Krylov solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();

  linearSolver->setPreconditioner(linearFickPreconditioner);

  nlinBVPOp->apply(rhsVec, solVec, resVec, 1.0, -1.0);
  double initialResidualNorm  = resVec->L2Norm();

  AMP::pout<<"Initial Residual Norm: "<<initialResidualNorm<<std::endl;

  nonlinearSolver->setZeroInitialGuess(false);

  nonlinearSolver->solve(rhsVec, solVec);

  nlinBVPOp->apply(rhsVec, solVec, resVec, 1.0, -1.0);

  double finalResidualNorm  = resVec->L2Norm();

  std::cout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

  solVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
  resVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // evaluate and register material coefficients for graphical output

  AMP::LinearAlgebra::Variable::shared_ptr fickCoeffVar(new AMP::Mesh::NodalScalarVariable("FickCoefficient"));
  AMP::LinearAlgebra::Variable::shared_ptr soretCoeffVar(new AMP::Mesh::NodalScalarVariable("SoretCoefficient"));
  AMP::LinearAlgebra::Vector::shared_ptr fickCoeffVec = meshAdapter->createVector(fickCoeffVar);
  AMP::LinearAlgebra::Vector::shared_ptr soretCoeffVec = meshAdapter->createVector(soretCoeffVar);
  manager->registerVectorAsData ( fickCoeffVec, "FickCoefficient" );
  manager->registerVectorAsData ( soretCoeffVec, "ThermalDiffusionCoefficient" );
  //boost::shared_ptr<AMP::Operator::DiffusionTransportModel> fickModel = fickOp->getTransportModel();
  //boost::shared_ptr<AMP::Operator::DiffusionTransportModel> soretModel = soretOp->getTransportModel();

  {
	  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator iterator = meshAdapter->beginOwnedNode();
	  size_t nnodes = fickCoeffVec->getLocalSize(), node;
	  std::vector<int> gids(nnodes);
	  std::vector<double> temp(nnodes), conc(nnodes), fickCoeff(nnodes), soretCoeff(nnodes), burn(nnodes);
	  for(node=0 ; iterator != meshAdapter->endOwnedNode(); iterator++ ) {
		  gids[node] = iterator->globalID();
		  node++;
	  }
	  AMP_INSIST(node==nnodes, "invalid count");
	  fickFrozen[AMP::Operator::Diffusion::TEMPERATURE]->getValuesByGlobalID(nnodes, &gids[0], &temp[0]);
	  solVec->getValuesByGlobalID(nnodes, &gids[0], &conc[0]);
	  // this is  used to plot the fick and soret coefficnets used.  commenting it out till someone finds out.
	  //fickModel->getTransport(fickCoeff, temp, conc, burn);  // This generates a compile error
	  //soretModel->getTransport(soretCoeff, temp, conc, burn);  // This generates a compiler error
	  //fickCoeffVec->setValuesByGlobalID(nnodes, &gids[0], &fickCoeff[0]);
	  //soretCoeffVec->setValuesByGlobalID(nnodes, &gids[0], &soretCoeff[0]);
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // write graphical output

#ifdef USE_SILO
  manager->writeFile<AMP::Mesh::SiloIO> ( exeName , 0 );
#endif

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // store result
  {
      AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator iterator = meshAdapter->beginOwnedNode();
      iterator = meshAdapter->beginOwnedNode();
      size_t numNodes = 0;
	  for(; iterator != meshAdapter->endOwnedNode(); iterator++ ) numNodes++;
	  results.resize(numNodes);

      iterator = meshAdapter->beginOwnedNode();
      size_t iNode=0;
	  for(; iterator != meshAdapter->endOwnedNode(); iterator++ ) {
		size_t gid = iterator->globalID();
		results[iNode] = solVec->getValueByGlobalID(gid);
		iNode++;
	  }
  }

  if(finalResidualNorm>1.0e-08) {
    ut->failure(exeName);
  } else {
    ut->passes("PetscSNES Solver successfully solves a nonlinear mechanics equation with Jacobian provided, FGMRES for Krylov");
  }
  ut->passes(exeName);
}


int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
	  std::vector<double>  results;
	  fickSoretTest(&ut, "testPetscSNESSolver-NonlinearFickSoret-cylinder-OxMSRZC09-1", results);
	  ut.passes("fick-soret quadratic external T");
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


