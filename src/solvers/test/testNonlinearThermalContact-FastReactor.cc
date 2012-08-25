#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <string>

#include "materials/Material.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"


#include "vectors/Variable.h"
#include "vectors/SimpleVector.h"
#include "vectors/Vector.h"

#include "ampmesh/SiloIO.h"
#include "ampmesh/MeshVariable.h"

#include "operators/MassLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/MassLinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/ElementOperationFactory.h"
#include "operators/ColumnBoundaryOperator.h"
#include "operators/FlowFrapconOperator.h"

#include "operators/map/Map3Dto1D.h"
#include "operators/map/Map1Dto3D.h"
#include "operators/map/MapOperatorParameters.h"
#include "operators/NeumannVectorCorrectionParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/DirichletMatrixCorrection.h"
#include "operators/DirichletVectorCorrection.h"
#include "operators/RobinMatrixCorrection.h"
#include "operators/RobinVectorCorrection.h"
#include "operators/NeumannVectorCorrection.h"
#include "operators/VolumeIntegralOperator.h"
#include "operators/NeutronicsRhs.h"

#include "../TrilinosMLSolver.h"
#include "../ColumnSolver.h"
#include "../PetscKrylovSolverParameters.h"
#include "../PetscKrylovSolver.h"
#include "../PetscSNESSolverParameters.h"
#include "../PetscSNESSolver.h"


void thermalContactTest(AMP::UnitTest *ut, std::string exeName )
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

//  AMP::AMPManager::startup();
//  AMP::Materials::initialize();

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::PIO::logAllNodes(log_file);

//  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
//  std::string mesh_file = input_db->getString("Mesh");

  AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( mgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter1 = manager->getMesh ( "pellet" );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter2 = manager->getMesh ( "clad" );

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  AMP::LinearAlgebra::Variable::shared_ptr TemperatureVar ( new AMP::Mesh::NodalScalarVariable ( "Temperature" ) );
  AMP::LinearAlgebra::Variable::shared_ptr inputVariable1 ( new AMP::Mesh::NodalScalarVariable ( "Temperature", meshAdapter1  ) );
  AMP::LinearAlgebra::Variable::shared_ptr inputVariable2 ( new AMP::Mesh::NodalScalarVariable ( "Temperature", meshAdapter2 ) );

  double intguess = input_db->getDoubleWithDefault("InitialGuess",400);

  AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvin = manager->createVector ( TemperatureVar );
  TemperatureInKelvin->setToScalar ( intguess );
  manager->registerVectorAsData ( TemperatureInKelvin , "Temperature" );


//-----------------------------------------------
//   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
//-----------------------------------------------

  AMP_INSIST( input_db->keyExists("NonlinearThermalOperator1"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel1;
  boost::shared_ptr<AMP::Database> nonlinearThermalDatabase1 = input_db->getDatabase("NonlinearThermalOperator1");
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator1 = boost::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1, nonlinearThermalDatabase1, thermalTransportModel1));

  // initialize the input variable
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator1 =
		  boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearThermalOperator1->getVolumeOperator());

  // initialize the output variable
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable1 = thermalVolumeOperator1->getOutputVariable();

  AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec1 = TemperatureInKelvin->subsetVectorForVariable ( inputVariable1 );
  AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec1       = meshAdapter1->createVector( outputVariable1 );
  AMP::LinearAlgebra::Vector::shared_ptr ResidualVec1            = meshAdapter1->createVector( outputVariable1 );

//-------------------------------------
//   CREATE THE LINEAR THERMAL OPERATOR 1 ----
//-------------------------------------

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel1;
  boost::shared_ptr<AMP::InputDatabase> bvpDatabase1 = boost::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase("LinearThermalOperator1"));
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator1 = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1, bvpDatabase1, transportModel1));

  //-------------------------------------
  //  CREATE THE NEUTRONICS SOURCE  //
  //-------------------------------------
  AMP_INSIST(input_db->keyExists("NeutronicsOperator"), "Key ''NeutronicsOperator'' is missing!");
  boost::shared_ptr<AMP::Database>  neutronicsOp_db = input_db->getDatabase("NeutronicsOperator");
  boost::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ));
  neutronicsParams->d_MeshAdapter = meshAdapter1;
  boost::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(new AMP::Operator::NeutronicsRhs( neutronicsParams ));

  AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar = neutronicsOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   SpecificPowerVec = meshAdapter1->createVector( SpecificPowerVar );

  neutronicsOperator->apply(nullVec, nullVec, SpecificPowerVec, 1., 0.);

  //----------------------------------------------------------
  //  Integrate Nuclear Rhs over Desnity * Volume //
  //----------------------------------------------------------

  AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
  boost::shared_ptr<AMP::Database> sourceDatabase = input_db->getDatabase("VolumeIntegralOperator");
  boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator = boost::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1, sourceDatabase, stransportModel));

  // Create the power (heat source) vector.
  AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   PowerInWattsVec = meshAdapter1->createVector( PowerInWattsVar );
  PowerInWattsVec->zero();

  // convert the vector of specific power to power for a given basis.
  sourceOperator->apply(nullVec, SpecificPowerVec, PowerInWattsVec, 1., 0.);

//--------------------------------------

  AMP_INSIST(input_db->keyExists("NonlinearSolver"),   "Key ''NonlinearSolver'' is missing!");

  boost::shared_ptr<AMP::Database> nonlinearSolver_db1 = input_db->getDatabase("NonlinearSolver"); 
  boost::shared_ptr<AMP::Database>  linearSolver_db1   = nonlinearSolver_db1->getDatabase("LinearSolver"); 

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the nonlinear solver
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams1(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db1));

  // change the next line to get the correct communicator out
  nonlinearSolverParams1->d_comm = globalComm;
  nonlinearSolverParams1->d_pOperator = nonlinearThermalOperator1;
  nonlinearSolverParams1->d_pInitialGuess = TemperatureInKelvinVec1;

  boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver1(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams1));

  //----------------------------------------------------------------------------------------------------------------------------------------------//

  boost::shared_ptr<AMP::Database> thermalPreconditioner_db1 = linearSolver_db1->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams1(new AMP::Solver::SolverStrategyParameters(thermalPreconditioner_db1));
  thermalPreconditionerParams1->d_pOperator = linearThermalOperator1;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalPreconditioner1(new AMP::Solver::TrilinosMLSolver(thermalPreconditionerParams1));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register the preconditioner with the Jacobian free Krylov solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver1 = nonlinearSolver1->getKrylovSolver();
  linearSolver1->setPreconditioner(linearThermalPreconditioner1);
  nonlinearThermalOperator1->apply(RightHandSideVec1, TemperatureInKelvinVec1, ResidualVec1, 1.0, -1.0);

//---------------------------------------------
//     CREATE THE CONTACT GAP OPERATOR 
//---------------------------------------------

  AMP_INSIST(input_db->keyExists("GapOperator"), "Key ''GapOperator'' is missing!");
  boost::shared_ptr<AMP::InputDatabase> gapDatabase       = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("GapOperator"));

  double heff = (gapDatabase)->getDouble("Convective_Coefficient");
  boost::shared_ptr<AMP::LinearAlgebra::Variable> gapVariable(new AMP::LinearAlgebra::Variable("Gap"));

//--------------------------------------------
//   CREATE THE LINEAR THERMAL OPERATOR 2 ----
//--------------------------------------------

  AMP_INSIST( input_db->keyExists("LinearThermalOperator2"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel2;
  boost::shared_ptr<AMP::Database> linearThermalDatabase2 = input_db->getDatabase("LinearThermalOperator2");
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator2 = boost::dynamic_pointer_cast<
    AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2, linearThermalDatabase2, thermalTransportModel2));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  boost::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> thermalVolumeOperator2 =
		  boost::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>(linearThermalOperator2->getVolumeOperator());

  // initialize the output variable
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable2 = thermalVolumeOperator2->getOutputVariable();

  AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec2 = TemperatureInKelvin->subsetVectorForVariable ( inputVariable2 );
  AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec2       = meshAdapter2->createVector( outputVariable2 );
  AMP::LinearAlgebra::Vector::shared_ptr ResidualVec2            = meshAdapter2->createVector( outputVariable2 );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams2 (new AMP::Solver::SolverStrategyParameters(linearSolver_db1));
  mlSolverParams2->d_pOperator = linearThermalOperator2;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver>         mlSolver2(new AMP::Solver::TrilinosMLSolver(mlSolverParams2));
  mlSolver2->setZeroInitialGuess(true);

//-------------------------------------

  AMP::LinearAlgebra::Vector::shared_ptr variableFluxVec1 = meshAdapter1->createVector( inputVariable1 );
  AMP::LinearAlgebra::Vector::shared_ptr scratchTempVec1  = meshAdapter1->createVector( inputVariable1 );
  variableFluxVec1->setToScalar(0.0);

  AMP::LinearAlgebra::Vector::shared_ptr variableFluxVec2 = meshAdapter2->createVector( inputVariable2 );
  AMP::LinearAlgebra::Vector::shared_ptr scratchTempVec2  = meshAdapter2->createVector( inputVariable2 );
  variableFluxVec2->setToScalar(0.0);

//-------------------------------------

  boost::shared_ptr<AMP::InputDatabase> map3dto1d_db1  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapPelletto1D"));
  boost::shared_ptr<AMP::Operator::MapOperatorParameters> map3dto1dParams1 (new AMP::Operator::MapOperatorParameters( map3dto1d_db1 ));
  map3dto1dParams1->d_MeshAdapter = meshAdapter1;
  boost::shared_ptr<AMP::Operator::Map3Dto1D> map1ToLowDim (new AMP::Operator::Map3Dto1D( map3dto1dParams1 ));

  boost::shared_ptr<AMP::InputDatabase> map1dto3d_db1  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("Map1DtoClad"));
  boost::shared_ptr<AMP::Operator::MapOperatorParameters> map1dto3dParams1 (new AMP::Operator::MapOperatorParameters( map1dto3d_db1 ));
  map1dto3dParams1->d_MapAdapter = meshAdapter2;
  //-------------------------------------
  // This is related to But # 1219 and 1210.
  //  -- It dies in compute_Z_locations of the constructor for mat1dto3d.  
  ut.passes("Everything up till constructing 1Dto3D passes.");
//if( AMP::AMP_MPI::getNodes() == 1 ) {
  //-------------------------------------
  boost::shared_ptr<AMP::Operator::Map1Dto3D> map1ToHighDim (new AMP::Operator::Map1Dto3D( map1dto3dParams1 ));

  map1ToLowDim->setZLocations( map1ToHighDim->getZLocations());

  boost::shared_ptr<AMP::InputDatabase> map3dto1d_db2  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapCladto1D"));
  boost::shared_ptr<AMP::Operator::MapOperatorParameters> map3dto1dParams2 (new AMP::Operator::MapOperatorParameters( map3dto1d_db2 ));
  map3dto1dParams2->d_MeshAdapter = meshAdapter2;
  boost::shared_ptr<AMP::Operator::Map3Dto1D> map2ToLowDim (new AMP::Operator::Map3Dto1D( map3dto1dParams2 ));

  boost::shared_ptr<AMP::InputDatabase> map1dto3d_db2  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("Map1DtoPellet"));
  boost::shared_ptr<AMP::Operator::MapOperatorParameters> map1dto3dParams2 (new AMP::Operator::MapOperatorParameters( map1dto3d_db2 ));
  map1dto3dParams2->d_MapAdapter = meshAdapter1;
  boost::shared_ptr<AMP::Operator::Map1Dto3D> map2ToHighDim (new AMP::Operator::Map1Dto3D( map1dto3dParams2 ));

  map2ToLowDim->setZLocations( map2ToHighDim->getZLocations());

///-------------------------------------
// From flow operator test  
//-------------------------------------

  boost::shared_ptr<AMP::InputDatabase> mapcladflow_db  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapCladtoFlow"));
  boost::shared_ptr<AMP::Operator::MapOperatorParameters> mapcladflowParams (new AMP::Operator::MapOperatorParameters( mapcladflow_db ));
  mapcladflowParams->d_MeshAdapter = meshAdapter2;
  boost::shared_ptr<AMP::Operator::Map3Dto1D> mapCladToFlow (new AMP::Operator::Map3Dto1D( mapcladflowParams ));

  boost::shared_ptr<AMP::InputDatabase> mapflowclad_db  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapFlowtoClad"));
  boost::shared_ptr<AMP::Operator::MapOperatorParameters> mapflowcladParams (new AMP::Operator::MapOperatorParameters( mapflowclad_db ));
  mapflowcladParams->d_MapAdapter = meshAdapter2;
  boost::shared_ptr<AMP::Operator::Map1Dto3D> mapFlowToClad (new AMP::Operator::Map1Dto3D( mapflowcladParams ));

  mapCladToFlow->setZLocations( mapFlowToClad->getZLocations());

//-------------------------------------
  size_t flowVecSize = mapFlowToClad->getNumZlocations(); 

//--------------------------------------
//     CREATE THE FLOW OPERATOR   ------
//--------------------------------------

  AMP_INSIST(input_db->keyExists("FlowFrapconOperator"), "Key ''FlowFrapconOperator'' is missing!");

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> flowtransportModel;
  boost::shared_ptr<AMP::InputDatabase> flowDatabase       = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("FlowFrapconOperator"));
  flowDatabase->putInteger("numPoints", flowVecSize);
  boost::shared_ptr<AMP::Operator::FlowFrapconOperator> flowOperator = boost::dynamic_pointer_cast<AMP::Operator::FlowFrapconOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2, flowDatabase, flowtransportModel));

  flowOperator->setZLocations( mapFlowToClad->getZLocations());

  AMP::LinearAlgebra::Variable::shared_ptr   inputVariable  =  flowOperator->getInputVariable() ;
  AMP::LinearAlgebra::Variable::shared_ptr   outputVariable =  flowOperator->getOutputVariable() ;

  double Cp, De, G, K, Re, Pr, hclad, dz, Tc, Tin;

  Cp   = (flowDatabase)->getDouble("Heat_Capacity");
  De   = (flowDatabase)->getDouble("Channel_Diameter");
  G    = (flowDatabase)->getDouble("Mass_Flux");
  K    = (flowDatabase)->getDouble("Conductivity");
  Re   = (flowDatabase)->getDouble("Reynolds");
  Pr   = (flowDatabase)->getDouble("Prandtl");
  Tin  = (flowDatabase)->getDouble("Temp_Inlet");

  hclad = (0.023*K/De)*pow(Re,0.8)*pow(Pr,0.4);

  dz = 2.5/flowVecSize ;

  AMP::LinearAlgebra::Vector::shared_ptr solVec  = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , inputVariable  );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec  = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , outputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr resVec  = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , outputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr workVec = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , inputVariable  );

  AMP::LinearAlgebra::Vector::shared_ptr vecLag  = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , outputVariable );

  resVec->setToScalar(350);
//-------------------------------------

  AMP::LinearAlgebra::Vector::shared_ptr robinRHSVec = meshAdapter2->createVector( thermalVolumeOperator2->getInputVariable() );

//------------------------------------------

  AMP::Operator::Operator::shared_ptr boundaryOp1;
  boundaryOp1 = nonlinearThermalOperator1->getBoundaryOperator(); 

  AMP::Operator::Operator::shared_ptr  robinBoundaryOp1;
//  robinBoundaryOp1 = (boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>(boundaryOp1) )->getBoundaryOperator(0);
  robinBoundaryOp1 = (boost::dynamic_pointer_cast<AMP::Operator::BoundaryOperator>(boundaryOp1) );

  boost::shared_ptr<AMP::InputDatabase> boundaryDatabase1      = boost::dynamic_pointer_cast<AMP::InputDatabase>( nonlinearThermalDatabase1->getDatabase("BoundaryOperator"));
//  boost::shared_ptr<AMP::InputDatabase> robinboundaryDatabase1 = boost::dynamic_pointer_cast<AMP::InputDatabase>( boundaryDatabase1->getDatabase("RobinVectorCorrection"));
  boost::shared_ptr<AMP::InputDatabase> robinboundaryDatabase1 = boost::dynamic_pointer_cast<AMP::InputDatabase>(boundaryDatabase1);

  robinboundaryDatabase1->putBool("constant_flux", false);
  robinboundaryDatabase1->putBool("skip_matrix_correction", true);
  boost::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters1 (new AMP::Operator::NeumannVectorCorrectionParameters( robinboundaryDatabase1 ) );

//------------------------------------------

  AMP::Operator::Operator::shared_ptr boundaryOp2;
  boundaryOp2 = linearThermalOperator2->getBoundaryOperator(); 
  AMP::Operator::Operator::shared_ptr  robinBoundaryOp2, robinBoundaryOp3;
  robinBoundaryOp2 = (boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>(boundaryOp2) )->getBoundaryOperator(0);
  robinBoundaryOp3 = (boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>(boundaryOp2) )->getBoundaryOperator(1);

  boost::shared_ptr<AMP::InputDatabase> boundaryDatabase2      = boost::dynamic_pointer_cast<AMP::InputDatabase>( linearThermalDatabase2->getDatabase("BoundaryOperator"));
  boost::shared_ptr<AMP::InputDatabase> robinboundaryDatabase2 = boost::dynamic_pointer_cast<AMP::InputDatabase>( boundaryDatabase2->getDatabase("RobinMatrixCorrection1"));
  boost::shared_ptr<AMP::InputDatabase> robinboundaryDatabase3 = boost::dynamic_pointer_cast<AMP::InputDatabase>( boundaryDatabase2->getDatabase("RobinMatrixCorrection2"));

  robinboundaryDatabase2->putBool("constant_flux", false);
  robinboundaryDatabase2->putBool("skip_matrix_correction", true);
  boost::shared_ptr<AMP::Operator::RobinMatrixCorrectionParameters> correctionParameters2 (new AMP::Operator::RobinMatrixCorrectionParameters( robinboundaryDatabase2 ) );

  robinboundaryDatabase3->putBool("constant_flux", false);
  robinboundaryDatabase3->putBool("skip_matrix_correction", true);
  boost::shared_ptr<AMP::Operator::RobinMatrixCorrectionParameters> correctionParameters3 (new AMP::Operator::RobinMatrixCorrectionParameters( robinboundaryDatabase3 ) );

//-------------------------------------

  size_t gapVecCladSize = map1ToHighDim->getNumZlocations(); 
  AMP::LinearAlgebra::Vector::shared_ptr gapVecClad = AMP::LinearAlgebra::SimpleVector::create( gapVecCladSize, gapVariable );

  size_t gapVecPelletSize = map2ToHighDim->getNumZlocations(); 
  AMP::LinearAlgebra::Vector::shared_ptr gapVecPellet = AMP::LinearAlgebra::SimpleVector::create( gapVecPelletSize, gapVariable );

  int cnt=0;
  AMP::LinearAlgebra::Vector::shared_ptr vecLag1           = meshAdapter1->createVector( outputVariable1 );
  vecLag1->copyVector(TemperatureInKelvinVec1);
  AMP::LinearAlgebra::Vector::shared_ptr vecLag2 	    = meshAdapter2->createVector( outputVariable2 );
  vecLag2->copyVector(TemperatureInKelvinVec2);

  bool testPassed = false;

  int maxIt = input_db->getIntegerWithDefault("max_iterations", 100);

  while ( cnt < maxIt )
  {
	  cnt++;

          //-----------------------------------------------
          RightHandSideVec1->zero();
          RightHandSideVec2->zero();

          RightHandSideVec1->copyVector(PowerInWattsVec);
          std::cout << "PowerInWattsVec norm  inside loop = " << RightHandSideVec1->L2Norm() <<"\n";

          map2ToLowDim->apply(nullVec,TemperatureInKelvinVec2,gapVecPellet ,1.0, 0.0);
          map2ToHighDim->apply(nullVec,gapVecPellet , scratchTempVec1,1.0, 0.0);

          scratchTempVec1->scale(heff);
          variableFluxVec1->copyVector(scratchTempVec1); 

          correctionParameters1->d_variableFlux = variableFluxVec1;
          robinBoundaryOp1->reset(correctionParameters1);

          std::cout<<"Variable flux1 norm inside loop : "<< variableFluxVec1->L2Norm() <<endl;

          nonlinearThermalOperator1->modifyRHSvector(RightHandSideVec1);
          nonlinearThermalOperator1->modifyInitialSolutionVector(TemperatureInKelvinVec1);
	  nonlinearSolver1->solve(RightHandSideVec1, TemperatureInKelvinVec1);
          nonlinearThermalOperator1->apply(RightHandSideVec1, TemperatureInKelvinVec1, ResidualVec1);

          std::cout<<"Norm of TemperatureInKelvinVec1: "<< TemperatureInKelvinVec1->L2Norm() <<endl;

          //------------------------------------------------------------

          mapCladToFlow->apply(nullVec,TemperatureInKelvinVec2,solVec,1.0,-1.0);
          int cnt2 = 0;
          while ( cnt2 < maxIt  )
          {
            cnt2++;
            flowOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);
            if ( abs(resVec->L2Norm()-vecLag->L2Norm()) < .000005*vecLag->L2Norm() )
              break;
            else
              std::cout << "for iteration cnt = "<<cnt <<" --> " << vecLag->L2Norm() << " " <<  resVec->L2Norm() << std::endl;

            std::cout<<"Intermediate Flow Solution " <<std::endl;
            for( int i=0; i<flowVecSize ; i++) {
              std::cout<<" @i : "<< i<<" is "<<resVec->getValueByLocalID(i) ;
            }
            std::cout<<std::endl;
            vecLag->copyVector(resVec);
          }

          mapFlowToClad->apply(nullVec, resVec, robinRHSVec,1.0,-1.0);

          robinRHSVec->scale(hclad);
          correctionParameters3->d_variableFlux = robinRHSVec; 
          robinBoundaryOp3->reset(correctionParameters3);
   
          //------------------------------------------------------------
          map1ToLowDim->apply(nullVec,TemperatureInKelvinVec1,gapVecClad ,1.0, 0.0);
    
          std::cout<<"Norm of solVec after map1toLowDim: "<< gapVecClad->L2Norm() <<endl;

          map1ToHighDim->apply(nullVec,gapVecClad , scratchTempVec2,1.0, 0.0);

          std::cout<<"Norm of scratch2: "<< scratchTempVec2->L2Norm() <<endl;

          scratchTempVec2->scale(heff);
          variableFluxVec2->copyVector(scratchTempVec2); 

          correctionParameters2->d_variableFlux = variableFluxVec2; 
          robinBoundaryOp2->reset(correctionParameters2);

          std::cout<<"Variable flux2 norm inside loop : "<< variableFluxVec2->L2Norm() <<endl;

          linearThermalOperator2->modifyRHSvector(RightHandSideVec2);
          linearThermalOperator2->apply(RightHandSideVec2, TemperatureInKelvinVec2, ResidualVec2);
          mlSolver2->solve(RightHandSideVec2, TemperatureInKelvinVec2);

          //------------------------------------------------------------

          std::cout<<"Residual Norm on Pellet after "<<cnt<<" iteration is : " << ResidualVec1->L2Norm() << std::endl;
          std::cout<<"Residual Norm on Clad after "<<cnt<<" iteration is : " << ResidualVec2->L2Norm() << std::endl;

          vecLag2->subtract(TemperatureInKelvinVec2,vecLag2);

//          if( nodes == 2 ) {
#ifdef USES_SILO
            manager->writeFile<AMP::Mesh::SiloIO> ( exeName , 0 );
#endif
//          }
          if (vecLag2->L2Norm() < 1.e-6 ) {
            testPassed = true;
            break;
          } else {
            std::cout << "for iteration cnt = "<<cnt <<" --> " << vecLag1->L2Norm() << " " <<  vecLag2->L2Norm() << std::endl;
          }
          std::cout<<std::endl;

    vecLag1->copyVector(TemperatureInKelvinVec1);
    vecLag2->copyVector(TemperatureInKelvinVec2);
  }

  //-------------------------------------

  if(testPassed)
  {
    ut.passes("Seggregated solve of Composite Operator using control loop of Nonlinear Thermal+Robin->Map->Gap->Map->Ninlinear Thermal+Robin .");
  }
  else
  {
    ITFAILS;
  }

//} else {
//  ut.expected_failure("parallel map3D-1D and map1D-3D fail in parallel, see bug #1219.");
//}
  input_db.reset();

  ut.passes(exeName);

//  AMP::AMPManager::shutdown();

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
    thermalContactTest(ut, "testNonlinearThermalContactPicard_FASTREACTOR");
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



