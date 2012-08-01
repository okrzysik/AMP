#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include "boost/shared_ptr.hpp"

#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"

#include "ampmesh/SiloIO.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorSelector.h"
#include "operators/ColumnOperator.h"
#include "operators/CoupledOperator.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/ElementOperationFactory.h"
#include "operators/OperatorBuilder.h"
#include "operators/SubchannelTwoEqNonlinearOperator.h"
#include "operators/SubchannelTwoEqLinearOperator.h"
#include "operators/SubchannelPhysicsModel.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/PetscSNESSolver.h"
#include "solvers/TrilinosMLSolver.h"

#include "SubchannelHelpers.h"

void SubchannelSolve(AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file = "output_" + exeName;
    std::string silo_name = exeName;
    AMP::PIO::logAllNodes(log_file);
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    globalComm.barrier();

    // Read the input file
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( input_file , input_db );
    input_db->printClassData(AMP::plog);

    // Get the Mesh database and create the mesh parameters
    boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(database));
    meshParams->setComm(globalComm);

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr  subchannelMesh = AMP::Mesh::Mesh::buildMesh(meshParams);

    AMP::Mesh::Mesh::shared_ptr xyFaceMesh;
    xyFaceMesh = subchannelMesh->Subset( getFaceIterator( subchannelMesh , 1 ) );

    int DofsPerFace =  2;
//    AMP::Discretization::DOFManager::shared_ptr faceDOFManager = AMP::Discretization::simpleDOFManager::create( xyFaceMesh, AMP::Mesh::Face, 1, DofsPerFace, true);
    AMP::Discretization::DOFManager::shared_ptr faceDOFManager = AMP::Discretization::simpleDOFManager::create( subchannelMesh, getFaceIterator( subchannelMesh , 1 ), getFaceIterator( subchannelMesh , 0 ), DofsPerFace);
    AMP::LinearAlgebra::Variable::shared_ptr flowVariable (new AMP::LinearAlgebra::Variable("Flow"));
    AMP::LinearAlgebra::Vector::shared_ptr flowSolVec = AMP::LinearAlgebra::createVector( faceDOFManager , flowVariable, true );

    boost::shared_ptr<AMP::Database> subchannelPhysics_db = input_db->getDatabase("SubchannelPhysicsModel");
    boost::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
    boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));

    boost::shared_ptr<AMP::Database> nonlinearOperator_db = input_db->getDatabase("SubchannelTwoEqNonlinearOperator");
    boost::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(new AMP::Operator::SubchannelOperatorParameters( nonlinearOperator_db ));
    subchannelOpParams->d_Mesh = xyFaceMesh ;
    subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel ;

    boost::shared_ptr<AMP::Operator::SubchannelTwoEqNonlinearOperator> twoEqNonlinearOperator (new AMP::Operator::SubchannelTwoEqNonlinearOperator(subchannelOpParams));

    boost::shared_ptr<AMP::Operator::SubchannelTwoEqLinearOperator> twoEqLinearOperator (new AMP::Operator::SubchannelTwoEqLinearOperator(subchannelOpParams));

    boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
    boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));
    nonlinearSolverParams->d_comm = globalComm;
    nonlinearSolverParams->d_pOperator = twoEqNonlinearOperator;
    nonlinearSolverParams->d_pInitialGuess = flowSolVec;
    boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver;
    linearSolver = nonlinearSolver->getKrylovSolver();

    boost::shared_ptr<AMP::Solver::TrilinosMLSolver> preconditioner;
    boost::shared_ptr<AMP::Database> preconditioner_db = nonlinearSolver_db->getDatabase("SubchannelPreconditioner");

    boost::shared_ptr<AMP::Solver::SolverStrategyParameters> preconditionerParams(new AMP::Solver::SolverStrategyParameters(preconditioner_db));
    preconditionerParams->d_pOperator = twoEqLinearOperator;
    preconditioner.reset(new AMP::Solver::TrilinosMLSolver(preconditionerParams));

    linearSolver->setPreconditioner(preconditioner);


#ifdef USE_SILO
    // Register the quantities to plot
    AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO );
    siloWriter->registerVector( flowSolVec, xyFaceMesh, AMP::Mesh::Face, "SubchannelFlowPressure" );
    siloWriter->writeFile( silo_name , 0 );
#endif


}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {
    SubchannelSolve(&ut, "testSubchannelSolve");
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


