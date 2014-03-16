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

#include "utils/Writer.h"
#include "vectors/Vector.h"
#include "operators/libmesh/MassLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/libmesh/MassLinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/libmesh/VolumeIntegralOperator.h"
#include "operators/subchannel/FlowFrapconOperator.h"
#include "operators/subchannel/FlowFrapconJacobian.h"
#include "operators/NeutronicsRhs.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/ElementOperationFactory.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/CoupledOperator.h"
#include "operators/CoupledOperatorParameters.h"

#include "operators/map/Map3Dto1D.h"
#include "operators/map/Map1Dto3D.h"
#include "operators/map/MapOperatorParameters.h"

#include "operators/boundary/ColumnBoundaryOperator.h"
#include "operators/boundary/NeumannVectorCorrectionParameters.h"
#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/libmesh/RobinMatrixCorrection.h"
#include "operators/boundary/RobinVectorCorrection.h"
#include "operators/boundary/libmesh/NeumannVectorCorrection.h"
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

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");

  AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( mgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "bar" );

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  double intguess = input_db->getDoubleWithDefault("InitialGuess",400);

  //-----------------------------------------------
  //   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
  //-----------------------------------------------
  AMP_INSIST( input_db->keyExists("NonlinearThermalOperator"), "key missing!" );
  boost::shared_ptr<AMP::Database> thermalNonlinear_db = input_db->getDatabase("NonlinearThermalOperator");
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> thermalNonlinearOperator = boost::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
											"NonlinearThermalOperator",
											input_db,
											thermalTransportModel));

  // initialize the input variable
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator =
    boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(thermalNonlinearOperator->getVolumeOperator());


  AMP::LinearAlgebra::Vector::shared_ptr globalSolVec = meshAdapter->createVector( thermalVolumeOperator->getOutputVariable() );
  AMP::LinearAlgebra::Vector::shared_ptr globalRhsVec = meshAdapter->createVector( thermalVolumeOperator->getOutputVariable() );
  AMP::LinearAlgebra::Vector::shared_ptr globalResVec = meshAdapter->createVector( thermalVolumeOperator->getOutputVariable() );

  globalSolVec->setToScalar(intguess); 

  manager->registerVectorAsData ( globalSolVec , "Temperature" );

  //-------------------------------------
  //   CREATE THE LINEAR THERMAL OPERATOR ----
  //-------------------------------------

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
  boost::shared_ptr<AMP::InputDatabase> bvp_db = boost::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase("LinearThermalOperator"));
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
  neutronicsParams->d_MeshAdapter = meshAdapter;
  boost::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(new AMP::Operator::NeutronicsRhs( neutronicsParams ));

  AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar = neutronicsOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   SpecificPowerVec = meshAdapter->createVector( SpecificPowerVar );

  neutronicsOperator->apply(nullVec, nullVec, SpecificPowerVec, 1., 0.);

  //----------------------------------------------------------
  //  Integrate Nuclear Rhs over Desnity * Volume //
  //----------------------------------------------------------

  AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
  boost::shared_ptr<AMP::Database> sourceDatabase = input_db->getDatabase("VolumeIntegralOperator");
  boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator = boost::dynamic_pointer_cast<
    AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
											  "VolumeIntegralOperator",
											  input_db,
											  stransportModel));

  // Create the power (heat source) vector.
  AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   PowerInWattsVec = meshAdapter->createVector( PowerInWattsVar );
  PowerInWattsVec->zero();

  // convert the vector of specific power to power for a given basis.
  sourceOperator->apply(nullVec, SpecificPowerVec, PowerInWattsVec, 1., 0.);

  //--------------------------------------
  AMP_INSIST(input_db->keyExists("NonlinearSolver"),   "Key ''NonlinearSolver'' is missing!");

  boost::shared_ptr<AMP::Database> nonlinearSolver_db1 = input_db->getDatabase("NonlinearSolver");
  boost::shared_ptr<AMP::Database>    linearSolver_db1 = nonlinearSolver_db1->getDatabase("LinearSolver");

  //-------------------------------------

  boost::shared_ptr<AMP::InputDatabase> map3dto1d_db  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("Map3Dto1D"));
  boost::shared_ptr<AMP::Operator::MapOperatorParameters> map3dto1dParams (new AMP::Operator::MapOperatorParameters( map3dto1d_db ));
  map3dto1dParams->d_MeshAdapter = meshAdapter;
  boost::shared_ptr<AMP::Operator::Map3Dto1D> mapToLowDim (new AMP::Operator::Map3Dto1D( map3dto1dParams ));

  boost::shared_ptr<AMP::InputDatabase> map1dto3d_db  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("Map1Dto3D"));
  boost::shared_ptr<AMP::Operator::MapOperatorParameters> map1dto3dParams (new AMP::Operator::MapOperatorParameters( map1dto3d_db ));
  map1dto3dParams->d_MapAdapter = meshAdapter;
  boost::shared_ptr<AMP::Operator::Map1Dto3D> mapToHighDim (new AMP::Operator::Map1Dto3D( map1dto3dParams ));

  mapToLowDim->setZLocations( mapToHighDim->getZLocations());

  //-------------------------------------
  //--------------------------------------
  //     CREATE THE FLOW OPERATOR   ------
  //--------------------------------------

  AMP_INSIST(input_db->keyExists("FlowFrapconOperator"), "Key ''FlowFrapconOperator'' is missing!");

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> flowtransportModel;
  boost::shared_ptr<AMP::InputDatabase> flowDatabase       = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("FlowFrapconOperator"));
  boost::shared_ptr<AMP::Operator::FlowFrapconOperator> flowOperator = boost::dynamic_pointer_cast<AMP::Operator::FlowFrapconOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																						      "FlowFrapconOperator",
																						      input_db,
																						      flowtransportModel));

  boost::shared_ptr<AMP::InputDatabase> flowJacDatabase    = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("FlowFrapconJacobian"));
  boost::shared_ptr<AMP::Operator::FlowFrapconJacobian> flowJacobian = boost::dynamic_pointer_cast<AMP::Operator::FlowFrapconJacobian>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																						      "FlowFrapconJacobian",
																						      input_db,
																						      flowtransportModel));

  flowOperator->setZLocations( mapToHighDim->getZLocations());
  flowJacobian->setZLocations( mapToHighDim->getZLocations());

  AMP::LinearAlgebra::Variable::shared_ptr   inputVariable  =  flowOperator->getInputVariable() ;
  AMP::LinearAlgebra::Variable::shared_ptr   outputVariable =  flowOperator->getOutputVariable() ;

  size_t flowVecSize = mapToHighDim->getNumZlocations(); 

  AMP::LinearAlgebra::Vector::shared_ptr flowSolVec  = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , inputVariable  );
  AMP::LinearAlgebra::Vector::shared_ptr flowRhsVec  = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , outputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr flowResVec  = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , outputVariable );

  AMP::LinearAlgebra::Vector::shared_ptr flowSolViewVec = AMP::LinearAlgebra::MultiVector::view( flowSolVec , globalComm );
  AMP::LinearAlgebra::Vector::shared_ptr flowRhsViewVec = AMP::LinearAlgebra::MultiVector::view( flowRhsVec , globalComm );
  AMP::LinearAlgebra::Vector::shared_ptr flowResViewVec = AMP::LinearAlgebra::MultiVector::view( flowResVec , globalComm );

  AMP::LinearAlgebra::Vector::shared_ptr cladVec = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , inputVariable );

  flowOperator->setVector(cladVec);
  flowJacobian->setVector(cladVec);

  double Cp, De, G, K, Re, Pr, heff, nP;

  Cp  = (flowDatabase)->getDouble("Heat_Capacity");
  De  = (flowDatabase)->getDouble("Channel_Diameter");
  G   = (flowDatabase)->getDouble("Mass_Flux");
  K   = (flowDatabase)->getDouble("Conductivity");
  Re  = (flowDatabase)->getDouble("Reynolds");
  Pr  = (flowDatabase)->getDouble("Prandtl");
  nP  = (flowDatabase)->getDouble("numpoints");
  NULL_USE(Cp);  NULL_USE(G);  NULL_USE(nP);

  heff = (0.023*K/De)*pow(Re,0.8)*pow(Pr,0.4);


  std::cout << "The flow Heff : "<< heff << std::endl;

//  flowSolVec->setToScalar(300);

  ut->passes("set up to the iterations passes.");

  AMP::LinearAlgebra::Vector::shared_ptr globalSolMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;
  globalSolMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( globalSolVec );
  globalSolMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( flowSolViewVec );

  AMP::LinearAlgebra::Vector::shared_ptr globalSolMultiVectorView = AMP::LinearAlgebra::MultiVector::view( globalSolMultiVector, globalComm );
  //---------------------------------------------------------------------------------------------------------------------//
  AMP::LinearAlgebra::Vector::shared_ptr globalRhsMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;
  globalRhsMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( globalRhsVec );
  globalRhsMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( flowRhsViewVec );

  AMP::LinearAlgebra::Vector::shared_ptr globalRhsMultiVectorView = AMP::LinearAlgebra::MultiVector::view( globalRhsMultiVector, globalComm );
  //---------------------------------------------------------------------------------------------------------------------//
  //---------------------------------------------------------------------------------------------------------------------//
  AMP::LinearAlgebra::Vector::shared_ptr globalResMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;
  globalResMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( globalResVec );
  globalResMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( flowResViewVec );

  AMP::LinearAlgebra::Vector::shared_ptr globalResMultiVectorView = AMP::LinearAlgebra::MultiVector::view( globalResMultiVector, globalComm );
  //---------------------------------------------------------------------------------------------------------------------//
  //---------------------------------------------------------------------------------------------------------------------//

  //-------------------------------------

  AMP::LinearAlgebra::Vector::shared_ptr robinRHSVec = meshAdapter->createVector( thermalNonlinearOperator->getOutputVariable() );

  mapToLowDim->setVector(cladVec);
  mapToHighDim->setVector( robinRHSVec );

  //-------------------------------------
  AMP::Operator::Operator::shared_ptr boundaryOp = thermalNonlinearOperator->getBoundaryOperator(); 

  boost::shared_ptr<AMP::Operator::RobinVectorCorrection> robinBoundaryOp = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>(   thermalNonlinearOperator->getBoundaryOperator() );
  boost::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters = boost::dynamic_pointer_cast<AMP::Operator::NeumannVectorCorrectionParameters>(robinBoundaryOp->getParameters()) ;

  robinBoundaryOp->setVariableFlux( robinRHSVec );

  //------------------------------------------------------------------
  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver");
  boost::shared_ptr<AMP::Database>    linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver");

  //-------------------------------------
  //Coupling Map to the Nonlinear Operators

  boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
  //-------------------------------------
  boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledNonlinearParams1(new AMP::Operator::CoupledOperatorParameters(tmp_db));
  coupledNonlinearParams1->d_MapOperator = mapToHighDim;
  coupledNonlinearParams1->d_BVPOperator = thermalNonlinearOperator;
  boost::shared_ptr<AMP::Operator::CoupledOperator>           coupledNonlinearOperator1(new AMP::Operator::CoupledOperator(coupledNonlinearParams1));

  boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledNonlinearParams2(new AMP::Operator::CoupledOperatorParameters(tmp_db));
  coupledNonlinearParams2->d_MapOperator = mapToLowDim;
  coupledNonlinearParams2->d_BVPOperator = flowOperator;
  boost::shared_ptr<AMP::Operator::CoupledOperator>           coupledNonlinearOperator2(new AMP::Operator::CoupledOperator(coupledNonlinearParams2));

  //Column of Coupled Operators
  boost::shared_ptr<AMP::Operator::OperatorParameters> columnNonlinearParams (new AMP::Operator::OperatorParameters(tmp_db));
  boost::shared_ptr<AMP::Operator::ColumnOperator>     columnNonlinearOperator(new AMP::Operator::ColumnOperator(columnNonlinearParams));
  columnNonlinearOperator->append(coupledNonlinearOperator1);
  columnNonlinearOperator->append(coupledNonlinearOperator2);

  boost::shared_ptr<AMP::Operator::OperatorParameters> columnLinearParams(new AMP::Operator::OperatorParameters(tmp_db));
  boost::shared_ptr<AMP::Operator::ColumnOperator>     coupledLinearOperator(new AMP::Operator::ColumnOperator(columnLinearParams));
  coupledLinearOperator->append(thermalLinearOperator);
  coupledLinearOperator->append(flowJacobian);

  //----------------------------------------------------------------//
  //-------------------------------------------------------------------------//
  // initialize the column preconditioner which is a diagonal block preconditioner
  boost::shared_ptr<AMP::Database>                 columnPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(new AMP::Solver::SolverStrategyParameters(columnPreconditioner_db));
  columnPreconditionerParams->d_pOperator = coupledLinearOperator;
  boost::shared_ptr<AMP::Solver::ColumnSolver>             columnPreconditioner(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

  boost::shared_ptr<AMP::Database>                 thermalPreconditioner_db = columnPreconditioner_db->getDatabase("pelletThermalPreconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams(new AMP::Solver::SolverStrategyParameters(thermalPreconditioner_db));
  thermalPreconditionerParams->d_pOperator = thermalLinearOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver>         thermalPreconditioner(new AMP::Solver::TrilinosMLSolver(thermalPreconditionerParams));

  //---------------
  boost::shared_ptr<AMP::Database> JacobianSolver_db = input_db->getDatabase("Flow1DSolver"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> flowSolverParams (new AMP::Solver::SolverStrategyParameters(JacobianSolver_db));
  flowSolverParams->d_pOperator = flowJacobian;
  boost::shared_ptr<AMP::Solver::Flow1DSolver>         flowJacobianSolver(new AMP::Solver::Flow1DSolver(flowSolverParams));

  columnPreconditioner->append(thermalPreconditioner);
  columnPreconditioner->append(flowJacobianSolver);

  //--------------------------------------------------------------------//
  // initialize the nonlinear solver
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm          = globalComm;
  nonlinearSolverParams->d_pOperator     = columnNonlinearOperator;
  nonlinearSolverParams->d_pInitialGuess = globalSolMultiVector ;
  boost::shared_ptr<AMP::Solver::PetscSNESSolver>  nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  // register the preconditioner with the Jacobian free Krylov solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();
  linearSolver->setPreconditioner(columnPreconditioner);

  //-------------------------------------
  nonlinearSolver->setZeroInitialGuess(false);


  globalRhsVec->zero();

  globalRhsVec->copyVector(PowerInWattsVec);
  std::cout << "PowerInWattsVec norm  inside loop = " << globalRhsVec->L2Norm() <<"\n";

  robinBoundaryOp->reset(correctionParameters);

  thermalNonlinearOperator->modifyRHSvector(globalRhsVec);
  thermalNonlinearOperator->modifyInitialSolutionVector(globalSolVec);

  columnNonlinearOperator->apply(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector, 1.0, -1.0);
  AMP::pout<<"Initial Residual Norm for Step is: "<<globalResMultiVector->L2Norm()<<std::endl;

  nonlinearSolver->solve(globalRhsMultiVectorView, globalSolMultiVectorView);

  columnNonlinearOperator->apply(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector, 1.0, -1.0);
  AMP::pout<<"Final   Residual Norm for Step is: "<<globalResMultiVector->L2Norm()<<std::endl;

  std::cout<<"Intermediate Flow Solution " <<std::endl;

  std::vector<double> expectedSolution(flowVecSize, 0);
  expectedSolution = (input_db->getDatabase("regression"))->getDoubleArray("expectedSolution");
  for(unsigned int i=0; i<flowVecSize ; i++) {
    if( !AMP::Utilities::approx_equal(expectedSolution[i], flowSolVec->getValueByLocalID(i), 1e-6)) {
      if(AMP::AMP_MPI::getRank() == 0) printf("solution: %.7e expected: %.7e \n", 
                                       flowSolVec->getValueByLocalID(i),expectedSolution[i] );
      ut->failure("solution is different for "+exeName);
    }
    // if(AMP::AMP_MPI::getRank() == 0) printf("%.7e, ",flowSolVec->getValueByLocalID(i) );
  }
  std::cout<<std::endl;

  //---------------------------------------------------------------------------

  if( AMP::AMP_MPI::getNodes() == 1 ) {
#ifdef USE_EXT_SILO
    manager->writeFile<AMP::Mesh::SiloIO> ( exeName , 0 );
#endif
  }

  if( globalResVec->L2Norm() < 10e-6 ) {
    ut->passes("Coupled solve of Composite Operator using control loop of Thermal+Robin->Map->Flow->Map .");
  } else {
    ut->failure("Coupled solve of Composite Operator using control loop of Thermal+Robin->Map->Flow->Map .");
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
        flowTest(&ut, "testThermalRobinFlow-1");
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



