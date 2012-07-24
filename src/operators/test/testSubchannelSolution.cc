#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "ampmesh/Mesh.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/SimpleVector.h"
#include "operators/SubchannelTwoEqNonlinearOperator.h"
#include "operators/OperatorBuilder.h"
#include "solvers/ColumnSolver.h"
#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/PetscSNESSolverParameters.h"
#include "solvers/PetscSNESSolver.h"
#include "solvers/TrilinosMLSolver.h"


void flowTest(AMP::UnitTest *ut, std::string exeName )
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;
  AMP::PIO::logAllNodes(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  // Read the input file
  boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
  AMP::InputManager::getManager()->parseInputFile ( input_file , input_db );

  // Get the Mesh database and create the mesh parameters
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(globalComm);

  // Create the meshes from the input database
  boost::shared_ptr<AMP::Mesh::Mesh> manager = AMP::Mesh::Mesh::buildMesh(meshParams);
  AMP::Mesh::Mesh::shared_ptr meshAdapter = manager->Subset( "bar" );

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  // get subchannel physics model
  boost::shared_ptr<AMP::Database> subchannelPhysics_db = input_db->getDatabase("SubchannelPhysicsModel");
  boost::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
  boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));

  // get nonlinear operator database
  boost::shared_ptr<AMP::Database> nonlinearOperator_db = input_db->getDatabase("SubchannelTwoEqNonlinearOperator");
  // create parameters
  boost::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(new AMP::Operator::SubchannelOperatorParameters( nonlinearOperator_db ));
  // put mesh into parameters
  subchannelOpParams->d_Mesh = meshAdapter;
  // put subchannel physics model into parameters
  subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;

  // create nonlinear operator
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
  boost::shared_ptr<AMP::Operator::SubchannelTwoEqNonlinearOperator> nonlinearOperator =
      boost::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqNonlinearOperator>(AMP::Operator::OperatorBuilder::createOperator(
      meshAdapter,"SubchannelTwoEqNonlinearOperator",input_db,elementModel ));
  
/*
  // create linear operator
  boost::shared_ptr<AMP::Operator::SubchannelTwoEqLinearOperator> linearOperator =
      boost::dynamic_pointer_cast<AMP::Operator::SubchannelTwoEqLinearOperator>(AMP::Operator::OperatorBuilder::createOperator(
      meshAdapter,"SubchannelTwoEqLinearOperator",input_db,elementModel ));
*/
  
  // pass creation test
  ut->passes(exeName+": creation");
  std::cout.flush();

  // get input and output variables
  AMP::LinearAlgebra::Variable::shared_ptr inputVariable  = nonlinearOperator->getInputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable = nonlinearOperator->getOutputVariable();

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr manufacturedVec = AMP::LinearAlgebra::SimpleVector::create( 5, inputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::SimpleVector::create( 5, inputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::SimpleVector::create( 5, outputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::SimpleVector::create( 5, outputVariable );

  // get exit pressure
  double Pout = nonlinearOperator_db->getDouble("Exit_Pressure");
  // set manufactured solution
  std::cout<<"Manufactured Solution:"<< std::endl;
  manufacturedVec->setValueByLocalID(0, 1000.0);
  manufacturedVec->setValueByLocalID(1, 15.3e6);
  manufacturedVec->setValueByLocalID(2, 15.2e6);
  manufacturedVec->setValueByLocalID(3, 15.1e6);
  manufacturedVec->setValueByLocalID(4, Pout);

  // get nonlinear solver database
  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
/*
  // get linear solver database
  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 
*/

  // put manufactured RHS into resVec
  nonlinearOperator->reset(subchannelOpParams);
  nonlinearOperator->apply(rhsVec, manufacturedVec, resVec, 1.0, 0.0);
/*
  linearOperator->reset(nonlinearOperator->getJacobianParameters(mv_view_solVec));
  linearOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);
*/
  
  // create nonlinear solver parameters
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

  // calculate initial guess
  // get inlet temperature
  double Tin = nonlinearOperator_db->getDouble("Inlet_Temperature");
  // compute inlet enthalpy
  std::map<std::string, boost::shared_ptr<std::vector<double> > > enthalpyArgMap;
  enthalpyArgMap.insert(std::make_pair("temperature",new std::vector<double>(1,Tin)));
  enthalpyArgMap.insert(std::make_pair("pressure",   new std::vector<double>(1,Pout)));
  std::vector<double> enthalpyResult(1);
  subchannelPhysicsModel->getProperty("Enthalpy",enthalpyResult,enthalpyArgMap); 
  double hin = enthalpyResult[0];
  // set initial guesses
  solVec->setValueByLocalID(0, hin);
  solVec->setValueByLocalID(1, Pout);
  solVec->setValueByLocalID(2, Pout);
  solVec->setValueByLocalID(3, Pout);
  solVec->setValueByLocalID(4, Pout);

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nonlinearOperator;
  nonlinearSolverParams->d_pInitialGuess = solVec;

/*
  // create linear solver parameters
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> linearSolverParams(new AMP::Solver::PetscSNESSolverParameters(linearSolver_db));

  // change the next line to get the correct communicator out
  linearSolverParams->d_comm = globalComm;
  linearSolverParams->d_pOperator = linearOperator;
  linearSolverParams->d_pInitialGuess = solVec;

  // create Jacobian solver
  boost::shared_ptr<AMP::Solver::PetscSNESSolver> JacobianSolver(new AMP::Solver::PetscSNESSolver(linearSolverParams));
*/

  // create nonlinear solver
  boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

/*
  // create linear solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();

  // create preconditioner
  boost::shared_ptr<AMP::Database> Preconditioner_db =  linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> PreconditionerParams(new AMP::Solver::SolverStrategyParameters(Preconditioner_db));
  PreconditionerParams->d_pOperator = linearOperator;
  boost::shared_ptr<AMP::Solver::Flow1DSolver> linearFlowPreconditioner(new AMP::Solver::Flow1DSolver(PreconditionerParams));
  // set preconditioner
  linearSolver->setPreconditioner(JacobianSolver);
*/

  // don't use zero initial guess
  nonlinearSolver->setZeroInitialGuess(false);

  // solve
  nonlinearSolver->solve(resVec, solVec);

  // print final solution
  for( int i=0; i<5; i++) {
     std::cout<<"Final_Solution["<<i<<"] = "<<solVec->getValueByLocalID(i);
     std::cout<<std::endl;
  }
  // print manufactured solution
  for( int i=0; i<5; i++) {
     std::cout<<"Manufactured_Solution["<<i<<"] = "<<manufacturedVec->getValueByLocalID(i);
     std::cout<<std::endl;
  }

  // print absolute error between manufactured solution and final solution
  AMP::LinearAlgebra::Vector::shared_ptr absErrorVec = AMP::LinearAlgebra::SimpleVector::create( 5, inputVariable );
  for( int i=0; i<5; i++) {
     absErrorVec->setValueByLocalID(i, manufacturedVec->getValueByLocalID(i) - solVec->getValueByLocalID(i));
     std::cout<<"Absolute_Error["<<i<<"] = "<<absErrorVec->getValueByLocalID(i);
     std::cout<<std::endl;
  }
  double absErrorNorm = absErrorVec->L2Norm();
  std::cout<<"L2 Norm of Absolute Error: "<<absErrorNorm<<std::endl;
  
  // print relative error between manufactured solution and final solution
  AMP::LinearAlgebra::Vector::shared_ptr relErrorVec = AMP::LinearAlgebra::SimpleVector::create( 5, inputVariable );
  for( int i=0; i<5; i++) {
     double relError;
     if (std::abs(manufacturedVec->getValueByLocalID(i)) < 1.0e-15){
        relError = 0.0;
     } else {
        relError = absErrorVec->getValueByLocalID(i)/manufacturedVec->getValueByLocalID(i);
     }
     relErrorVec->setValueByLocalID(i, relError);
     std::cout<<"Relative_Error["<<i<<"] = "<<relErrorVec->getValueByLocalID(i);
     std::cout<<std::endl;
  }
  double relErrorNorm = relErrorVec->L2Norm();
  std::cout<<"L2 Norm of Relative Error: "<<relErrorNorm<<std::endl;

  // check that norm of relative error is less than tolerance
  if(relErrorNorm > 1.0e-7){
     ut->failure(exeName+": manufactured solution test");
  } else {
     ut->passes(exeName+": manufactured solution test");
  }

  input_db.reset();

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
        flowTest(&ut, "testSubchannelSolution");
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

