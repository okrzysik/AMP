#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "materials/Material.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "vectors/Variable.h"

#include "ampmesh/Mesh.h"
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "ampmesh/SiloIO.h"
#include "vectors/Vector.h"

#include "operators/MassLinearElement.h"
#include "operators/MassLinearFEOperator.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionTransportModel.h"

#include "operators/VolumeIntegralOperator.h"
#include "operators/subchannel/FlowFrapconOperator.h"
#include "operators/subchannel/FlowFrapconJacobian.h"
#include "operators/NeutronicsRhs.h"

#include "operators/ElementPhysicsModelFactory.h"
#include "operators/ElementOperationFactory.h"

#include "operators/boundary/ColumnBoundaryOperator.h"
#include "operators/boundary/NeumannVectorCorrectionParameters.h"
#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/RobinMatrixCorrection.h"
#include "operators/boundary/RobinVectorCorrection.h"
#include "operators/boundary/NeumannVectorCorrection.h"

#include "operators/NonlinearBVPOperator.h"
#include "operators/CoupledOperator.h"
#include "operators/CoupledOperatorParameters.h"

#include "operators/map/Map3Dto1D.h"
#include "operators/map/Map1Dto3D.h"
#include "operators/map/MapOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "../TrilinosMLSolver.h"
#include "../ColumnSolver.h"
#include "../PetscKrylovSolverParameters.h"
#include "../PetscKrylovSolver.h"
#include "../PetscSNESSolverParameters.h"
#include "../PetscSNESSolver.h"
#include "solvers/Flow1DSolver.h"


void flowTest(AMP::UnitTest *ut, std::string exeName )
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::PIO::logAllNodes(log_file);
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);

//--------------------------------------------------
//   Create the Mesh.
//--------------------------------------------------
    AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
    boost::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
    boost::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
    mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
    boost::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh(mgrParams);
//--------------------------------------------------

//--------------------------------------------------
// Create a DOF manager for a nodal vector 
//--------------------------------------------------
  int DOFsPerNode = 1;
  int DOFsPerElement = 8;
  int nodalGhostWidth = 1;
  int gaussPointGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr nodalDofMap      = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth,      DOFsPerNode,    split);
  AMP::Discretization::DOFManager::shared_ptr gaussPointDofMap = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Volume, gaussPointGhostWidth, DOFsPerElement, split);
//--------------------------------------------------

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  double intguess = input_db->getDoubleWithDefault("InitialGuess",400);

  //-----------------------------------------------
  //   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
  //-----------------------------------------------
  AMP_INSIST( input_db->keyExists("NonlinearThermalOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> thermalNonlinearOperator = boost::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
											"NonlinearThermalOperator",
											input_db,
											thermalTransportModel));

  // initialize the input variable
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator =
    boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(thermalNonlinearOperator->getVolumeOperator());


  AMP::LinearAlgebra::Vector::shared_ptr globalSolVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVolumeOperator->getOutputVariable() );
  AMP::LinearAlgebra::Vector::shared_ptr globalRhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVolumeOperator->getOutputVariable() );
  AMP::LinearAlgebra::Vector::shared_ptr globalResVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVolumeOperator->getOutputVariable() );

  globalSolVec->setToScalar(intguess); 

  //-------------------------------------
  //   CREATE THE LINEAR THERMAL OPERATOR ----
  //-------------------------------------

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> thermalLinearOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																							   "LinearThermalOperator",
																							   input_db,
																							   transportModel));

  //-------------------------------------
  //  CREATE THE NEUTRONICS SOURCE  //
  //-------------------------------------
  AMP_INSIST(input_db->keyExists("NeutronicsOperator"), "Key ''NeutronicsOperator'' is missing!");
  boost::shared_ptr<AMP::Database>  neutronicsOp_db = input_db->getDatabase("NeutronicsOperator");
  boost::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ));
  boost::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(new AMP::Operator::NeutronicsRhs( neutronicsParams ));

  AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar = neutronicsOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   SpecificPowerVec = AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

  neutronicsOperator->apply(nullVec, nullVec, SpecificPowerVec, 1., 0.);

  //----------------------------------------------------------
  //  Integrate Nuclear Rhs over Desnity * Volume //
  //----------------------------------------------------------

  AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
  boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator = boost::dynamic_pointer_cast<
    AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
											  "VolumeIntegralOperator",
											  input_db,
											  stransportModel));

  // Create the power (heat source) vector.
  AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   PowerInWattsVec = AMP::LinearAlgebra::createVector( nodalDofMap, PowerInWattsVar );
  PowerInWattsVec->zero();

  // convert the vector of specific power to power for a given basis.
  sourceOperator->apply(nullVec, SpecificPowerVec, PowerInWattsVec, 1., 0.);

  //--------------------------------------
  AMP_INSIST(input_db->keyExists("NonlinearSolver"),   "Key ''NonlinearSolver'' is missing!");

  //boost::shared_ptr<AMP::Database> nonlinearSolver_db1 = input_db->getDatabase("NonlinearSolver");
  //boost::shared_ptr<AMP::Database>    linearSolver_db1 = nonlinearSolver_db1->getDatabase("LinearSolver");

  ut->passes("set up to the iterations passes.");
  //-------------------------------------

  AMP::LinearAlgebra::Vector::shared_ptr globalSolMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;
  globalSolMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( globalSolVec );

  //AMP::LinearAlgebra::Vector::shared_ptr globalSolMultiVectorView = AMP::LinearAlgebra::MultiVector::view( globalSolMultiVector, globalComm );
  //---------------------------------------------------------------------------------------------------------------------//
  AMP::LinearAlgebra::Vector::shared_ptr globalRhsMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;
  globalRhsMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( globalRhsVec );

  //AMP::LinearAlgebra::Vector::shared_ptr globalRhsMultiVectorView = AMP::LinearAlgebra::MultiVector::view( globalRhsMultiVector, globalComm );
  //---------------------------------------------------------------------------------------------------------------------//
  AMP::LinearAlgebra::Vector::shared_ptr globalResMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;
  globalResMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( globalResVec );

  //AMP::LinearAlgebra::Vector::shared_ptr globalResMultiVectorView = AMP::LinearAlgebra::MultiVector::view( globalResMultiVector, globalComm );
  //---------------------------------------------------------------------------------------------------------------------//


  //AMP::LinearAlgebra::Vector::shared_ptr robinRHSVec = meshAdapter->createVector( thermalNonlinearOperator->getOutputVariable() );

  //-------------------------------------
  AMP::Operator::Operator::shared_ptr boundaryOp = thermalNonlinearOperator->getBoundaryOperator(); 

  //  boost::shared_ptr<AMP::Operator::RobinVectorCorrection> robinBoundaryOp = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>(   thermalNonlinearOperator->getBoundaryOperator() );
  //  boost::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters = boost::dynamic_pointer_cast<AMP::Operator::NeumannVectorCorrectionParameters>(robinBoundaryOp->getParameters()) ;

  //  robinBoundaryOp->setVariableFlux( robinRHSVec );

  //------------------------------------------------------------------
  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver");
  boost::shared_ptr<AMP::Database>    linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver");

  //----------------------------------------------------------------//
  // initialize the nonlinear solver
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm          = globalComm;
  nonlinearSolverParams->d_pOperator     = thermalNonlinearOperator;
  nonlinearSolverParams->d_pInitialGuess = globalSolVec ;
  boost::shared_ptr<AMP::Solver::PetscSNESSolver>  nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));
  //-------------------------------------------------------------------------//
  // initialize the column preconditioner which is a diagonal block preconditioner
  boost::shared_ptr<AMP::Database>                 columnPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");

  boost::shared_ptr<AMP::Database>                 thermalPreconditioner_db = columnPreconditioner_db->getDatabase("pelletThermalPreconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams(new AMP::Solver::SolverStrategyParameters(thermalPreconditioner_db));
  thermalPreconditionerParams->d_pOperator = thermalLinearOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver>         thermalPreconditioner(new AMP::Solver::TrilinosMLSolver(thermalPreconditionerParams));

  //--------------------------------------------------------------------//
  // register the preconditioner with the Jacobian free Krylov solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();
  linearSolver->setPreconditioner(thermalPreconditioner);

  //-------------------------------------
  nonlinearSolver->setZeroInitialGuess(false);


  globalRhsVec->zero();

  globalRhsVec->copyVector(PowerInWattsVec);
  std::cout << "PowerInWattsVec norm  inside loop = " << globalRhsVec->L2Norm() <<"\n";
  double expectedVal = 0.175811;
  if( !AMP::Utilities::approx_equal( expectedVal, globalRhsVec->L2Norm(), 1e-5) ) { 
        ut->failure("the PowerInWattsVec norm has changed."); }

  //    robinBoundaryOp->reset(correctionParameters);

  thermalNonlinearOperator->modifyRHSvector(globalRhsVec);
  thermalNonlinearOperator->modifyInitialSolutionVector(globalSolVec);

  thermalNonlinearOperator->apply(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector, 1.0, -1.0);
  AMP::pout<<"Initial Residual Norm for Step is: "<<globalResVec->L2Norm()<<std::endl;
  expectedVal = 4.84311;
  if( !AMP::Utilities::approx_equal( expectedVal, globalResVec->L2Norm(), 1e-5) ) { 
        ut->failure("the Initial Residual Norm has changed."); }

  std::cout << " RHS Vec L2 Norm "<< globalRhsVec->L2Norm()<<std::endl;
  nonlinearSolver->solve(globalRhsVec, globalSolVec);

  std::cout<<"Final Solution Norm: "<<globalSolVec->L2Norm()<<std::endl;
  expectedVal = 51541;
  if( !AMP::Utilities::approx_equal( expectedVal, globalSolVec->L2Norm(), 1e-5) ) {
        ut->failure("the Final Solution Norm has changed."); }

  thermalNonlinearOperator->apply(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector, 1.0, -1.0);
  AMP::pout<<"Final   Residual Norm for Step is: "<<globalResVec->L2Norm()<<std::endl;
  expectedVal = 1.-10;
  if( !AMP::Utilities::approx_equal( expectedVal, globalResVec->L2Norm(), 10) ) {
        ut->failure("the Final Residual Norm has changed."); }

  //---------------------------------------------------------------------------

#ifdef USE_SILO
     AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
     siloWriter->registerMesh( meshAdapter );

     siloWriter->registerVector( globalSolVec, meshAdapter, AMP::Mesh::Vertex, "Temperature" );
 
     siloWriter->writeFile( input_file , 0 );
#endif

  if( globalResVec->L2Norm() < 10e-6 ) {
    ut->passes("Seggregated solve of Composite Operator using control loop of Thermal+Robin->Map->Flow->Map .");
  } else {
    ut->failure("Seggregated solve of Composite Operator using control loop of Thermal+Robin->Map->Flow->Map .");
  }


  //-------------------------------------
  // The 3D-to-1D map is not working in parallel.
  //   -- See Bug 1219 and 1209.
  //} else {
  //  ut.expected_failure("parallel map3D-1D and map1D-3D fail in parallel, see bug #1219.");
  //}
input_db.reset();

ut->passes(exeName);

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
        flowTest(&ut, "testThermalRobinFlow-2");
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



