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

#include "operators/map/MapSurface.h"
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
#include "operators/CoupledOperator.h"
#include "operators/CoupledOperatorParameters.h"

#include "solvers/PetscKrylovSolver.h"
#include "solvers/TrilinosMLSolver.h"
#include "solvers/ColumnSolver.h"
#include "time_integrators/ColumnTimeOperator.h"
#include "time_integrators/sundials/IDATimeIntegrator.h"
#include "time_integrators/sundials/IDATimeOperator.h"

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

void thermalContactTest(AMP::UnitTest *ut, std::string exeName )
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

//  AMP::AMPManager::startup();
//  

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  AMP::PIO::logAllNodes(log_file);

//  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
//  std::string mesh_file = input_db->getString("Mesh");

  AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( mgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter1 = manager->getMesh ( "pellet" );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter2 = manager->getMesh ( "clad" );

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  AMP::LinearAlgebra::Variable::shared_ptr TemperatureVar ( new AMP::NodalScalarVariable ( "Temperature" ) );
  AMP::LinearAlgebra::Variable::shared_ptr inputVariable1 ( new AMP::NodalScalarVariable ( "Temperature", meshAdapter1 ) );
  AMP::LinearAlgebra::Variable::shared_ptr inputVariable2 ( new AMP::NodalScalarVariable ( "Temperature", meshAdapter2 ) );

  AMP::LinearAlgebra::Variable::shared_ptr outputVariable1 ( new AMP::NodalScalarVariable ( "Temperature", meshAdapter1 ) );
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable2 ( new AMP::NodalScalarVariable ( "Temperature", meshAdapter2 ) );

  double initialGuess = input_db->getDoubleWithDefault("InitialGuess",400);

  AMP::LinearAlgebra::Vector::shared_ptr temperature = manager->createVector ( TemperatureVar );
  AMP::LinearAlgebra::Vector::shared_ptr temperatureICPrime = manager->createVector ( TemperatureVar );
  AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec    = manager->createVector ( TemperatureVar );    
  AMP::LinearAlgebra::Vector::shared_ptr ResidualVec         = manager->createVector ( TemperatureVar );    
  AMP::LinearAlgebra::Vector::shared_ptr scratchVec           = manager->createVector ( TemperatureVar );    
  AMP::LinearAlgebra::Vector::shared_ptr temperatureVec1 = temperature->subsetVectorForVariable ( inputVariable1  );
  AMP::LinearAlgebra::Vector::shared_ptr temperatureICPrime1 = temperatureICPrime->subsetVectorForVariable ( inputVariable1  );
  AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec1       = RightHandSideVec->subsetVectorForVariable    ( outputVariable1 );
  AMP::LinearAlgebra::Vector::shared_ptr ResidualVec1            = ResidualVec->subsetVectorForVariable         ( outputVariable1 );
  AMP::LinearAlgebra::Vector::shared_ptr scratchVec1           = scratchVec->subsetVectorForVariable         ( outputVariable1 );
  AMP::LinearAlgebra::Vector::shared_ptr temperatureVec2 = temperature->subsetVectorForVariable ( inputVariable2  );
  AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec2       = RightHandSideVec->subsetVectorForVariable    ( outputVariable2 );
  AMP::LinearAlgebra::Vector::shared_ptr ResidualVec2            = ResidualVec->subsetVectorForVariable         ( outputVariable2 );


  temperature->setToScalar ( initialGuess );
  temperatureICPrime->zero();
  
  manager->registerVectorAsData ( temperature , "Temperature" );
  manager->registerVectorAsData ( ResidualVec         , "Residual" );

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
  //  Integrate Nuclear Rhs over Density * Volume //
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

  //  PowerInWattsVec->scale(1.0e-03);
  //copy the power to pellet RHS vector
  RightHandSideVec1->copyVector(PowerInWattsVec);

  cout<<"L2 Norm of the PowerInWattsVec and RightHSVec1"<<PowerInWattsVec->L2Norm()<<" "<<RightHandSideVec1->L2Norm()<<endl;


  //-----------------------------------------------
  //   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
  //-----------------------------------------------
  
  AMP_INSIST( input_db->keyExists("NonlinearThermalOperator1"), "key missing!" );
  
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel1;
  boost::shared_ptr<AMP::Database> nonlinearThermalDatabase1 = input_db->getDatabase("NonlinearThermalOperator1");
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator1 = boost::dynamic_pointer_cast<
  AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1, nonlinearThermalDatabase1, thermalTransportModel1));
  
  //-------------------------------------
  //   CREATE THE LINEAR THERMAL OPERATOR 1 ----
  //-------------------------------------
  
  boost::shared_ptr<AMP::InputDatabase> bvpDatabase1 = boost::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase("LinearThermalOperator1"));
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator1 = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1, bvpDatabase1, thermalTransportModel1));
  
  //--------------------------------------------
  //   CREATE THE NONLINEAR THERMAL OPERATOR 2 ----
  //--------------------------------------------
  
  AMP_INSIST( input_db->keyExists("NonlinearThermalOperator2"), "key missing!" );
  
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel2;
  boost::shared_ptr<AMP::Database> nonlinearThermalDatabase2 = input_db->getDatabase("NonlinearThermalOperator2");
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator2 = boost::dynamic_pointer_cast<
  AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2, nonlinearThermalDatabase2, thermalTransportModel2));
  
  //--------------------------------------------
  //   CREATE THE LINEAR THERMAL OPERATOR 2 ----
  //--------------------------------------------
  
  boost::shared_ptr<AMP::InputDatabase> bvpDatabase2 = boost::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase("LinearThermalOperator2"));
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator2 = boost::dynamic_pointer_cast<
  AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2, bvpDatabase2, thermalTransportModel2));
  
  //-------------------------------------
  // initialize the mapping operators
  
  boost::shared_ptr<AMP::InputDatabase> mapcladtopellet_db = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapCladtoPellet")); 
  boost::shared_ptr<AMP::Operator::MapSurface> mapcladtopellet = boost::dynamic_pointer_cast<AMP::Operator::MapSurface>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2, meshAdapter1, mapcladtopellet_db ));
  
  boost::shared_ptr<AMP::InputDatabase> mappellettoclad_db = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapPellettoClad")); 
  boost::shared_ptr<AMP::Operator::MapSurface> mappellettoclad = boost::dynamic_pointer_cast<AMP::Operator::MapSurface>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1, meshAdapter2, mappellettoclad_db ));
  
  //------------------------------------------
  // initialize the input variable
  
  AMP::Operator::Operator::shared_ptr  robinBoundaryOp1 = nonlinearThermalOperator1->getBoundaryOperator(); 

  boost::shared_ptr<AMP::InputDatabase> robinboundaryDatabase1 = boost::dynamic_pointer_cast<AMP::InputDatabase>( nonlinearThermalDatabase1->getDatabase("BoundaryOperator"));

  robinboundaryDatabase1->putBool("constant_flux", false);
  robinboundaryDatabase1->putBool("skip_matrix_correction", true);
  boost::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters1 (new AMP::Operator::NeumannVectorCorrectionParameters( robinboundaryDatabase1 ) );

//------------------------------------------

  AMP::Operator::Operator::shared_ptr boundaryOp2;
  boundaryOp2 = nonlinearThermalOperator2->getBoundaryOperator(); 

  AMP::Operator::Operator::shared_ptr  robinBoundaryOp2;
  robinBoundaryOp2 = (boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>(boundaryOp2) )->getBoundaryOperator(0);

  boost::shared_ptr<AMP::InputDatabase> boundaryDatabase2      = boost::dynamic_pointer_cast<AMP::InputDatabase>( nonlinearThermalDatabase2->getDatabase("BoundaryOperator"));
  boost::shared_ptr<AMP::InputDatabase> robinboundaryDatabase2 = boost::dynamic_pointer_cast<AMP::InputDatabase>( boundaryDatabase2->getDatabase("RobinVectorCorrection"));

  robinboundaryDatabase2->putBool("constant_flux", false);
  robinboundaryDatabase2->putBool("skip_matrix_correction", true);
  boost::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters2 (new AMP::Operator::NeumannVectorCorrectionParameters( robinboundaryDatabase2 ) );

//-------------------------------------
  //Coupling Map to the Nonlinear Operators
  boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
  boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledNonlinearPelletParams(new AMP::Operator::CoupledOperatorParameters(tmp_db));

  coupledNonlinearPelletParams->d_MapOperator = mapcladtopellet;
  coupledNonlinearPelletParams->d_BVPOperator = nonlinearThermalOperator1;
  boost::shared_ptr<AMP::Operator::CoupledOperator> coupledNonlinearPellet(new AMP::Operator::CoupledOperator(coupledNonlinearPelletParams));
  //-------------------------------------
  boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledNonlinearCladParams(new AMP::Operator::CoupledOperatorParameters(tmp_db));
  coupledNonlinearCladParams->d_MapOperator = mappellettoclad ;
  coupledNonlinearCladParams->d_BVPOperator = nonlinearThermalOperator2;
  boost::shared_ptr<AMP::Operator::CoupledOperator> coupledNonlinearClad(new AMP::Operator::CoupledOperator(coupledNonlinearCladParams));

  //-------------------------------------
  //Column of nonlinear coupled operators
  boost::shared_ptr<AMP::Operator::OperatorParameters> nonlinearParams (new AMP::Operator::OperatorParameters(tmp_db));
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearCoupledOperator(new AMP::Operator::ColumnOperator(nonlinearParams));
  nonlinearCoupledOperator->append(coupledNonlinearPellet);
  nonlinearCoupledOperator->append(coupledNonlinearClad);

  //-------------------------------------
  //Column of linear operators
  boost::shared_ptr<AMP::Operator::OperatorParameters> linearParams(new AMP::Operator::OperatorParameters(tmp_db));
  boost::shared_ptr<AMP::Operator::ColumnOperator> linearCoupledOperator(new AMP::Operator::ColumnOperator(linearParams));
  linearCoupledOperator->append(linearThermalOperator1);
  linearCoupledOperator->append(linearThermalOperator2);

  // ---------------------------------------------------------------------------------------
  // create a mass linear BVP operator
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> massThermalModel1;
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> massThermalModel2;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> massThermalOp1;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> massThermalOp2;
  boost::shared_ptr<AMP::InputDatabase> massThermalOp_db1 =  boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MassThermalOperator1"));
  boost::shared_ptr<AMP::InputDatabase> massThermalOp_db2 =  boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MassThermalOperator2"));
  massThermalOp1 = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1, massThermalOp_db1,massThermalModel1));
  massThermalOp2 = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2, massThermalOp_db2,massThermalModel2));
  
  // ---------------------------------------------------------------------------------------
  // create a column mass operator object for use in the nonlinear and linear problem definition
  boost::shared_ptr<AMP::Operator::OperatorParameters> massParams;
  boost::shared_ptr<AMP::Operator::ColumnOperator> columnMassOperator(new AMP::Operator::ColumnOperator(massParams));
  columnMassOperator->append(massThermalOp1);
  columnMassOperator->append(massThermalOp2);

  //--------------------------------------------------------------------//

  AMP::LinearAlgebra::Variable::shared_ptr TemperatureMapVar ( new AMP::NodalScalarVariable ( "TemperatureMap" ) );
  AMP::LinearAlgebra::Variable::shared_ptr inputMapVariable1 ( new AMP::NodalScalarVariable ( "TemperatureMap", meshAdapter1 ) );
  AMP::LinearAlgebra::Variable::shared_ptr inputMapVariable2 ( new AMP::NodalScalarVariable ( "TemperatureMap", meshAdapter2 ) );

  AMP::LinearAlgebra::Vector::shared_ptr TemperatureMapvector  = manager->createVector ( TemperatureMapVar );
  AMP::LinearAlgebra::Vector::shared_ptr scratchCladToPellet  = TemperatureMapvector->subsetVectorForVariable ( inputMapVariable1  );
  AMP::LinearAlgebra::Vector::shared_ptr scratchPelletToClad  = TemperatureMapvector->subsetVectorForVariable ( inputMapVariable2  );

  correctionParameters1->d_variableFlux = scratchCladToPellet ;
  correctionParameters2->d_variableFlux = scratchPelletToClad ;

  robinBoundaryOp1->reset(correctionParameters1);
  robinBoundaryOp2->reset(correctionParameters2);

  mapcladtopellet->setVector(scratchCladToPellet ); 
  mappellettoclad->setVector(scratchPelletToClad ); 

  //-------------------------------------
  //Applying boundary conditions to the nonlinear BVP Operator

  nonlinearThermalOperator1->modifyRHSvector(RightHandSideVec1);
  nonlinearThermalOperator1->modifyInitialSolutionVector(temperatureVec1);

  nonlinearThermalOperator2->modifyRHSvector(RightHandSideVec2);
  nonlinearThermalOperator2->modifyInitialSolutionVector(temperatureVec2);

  //-------------------------------------------------------------------------------------------
  // solve for the time derivative at the initial time
  // form the rhs for the solve
  nonlinearCoupledOperator->apply(RightHandSideVec, temperature, scratchVec, -1.0, 1.0);
#if 0

  boost::shared_ptr<AMP::Database> icLinearSolver_db = input_db->getDatabase("ICLinearSolver"); 
  // ---- first initialize the preconditioner
  boost::shared_ptr<AMP::Database> icpcSolver_db = icLinearSolver_db->getDatabase("Preconditioner"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> icColumnPreconditionerParams(new AMP::Solver::SolverStrategyParameters(icpcSolver_db));
  icColumnPreconditionerParams->d_pOperator = columnMassOperator;
  boost::shared_ptr<AMP::Solver::ColumnSolver> icColumnPreconditioner(new AMP::Solver::ColumnSolver(icColumnPreconditionerParams));

  boost::shared_ptr<AMP::Database> icPelletPreconditioner_db = icpcSolver_db->getDatabase("pelletPreconditioner"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> icPelletPreconditionerParams(new AMP::Solver::SolverStrategyParameters(icPelletPreconditioner_db));
  icPelletPreconditionerParams->d_pOperator = columnMassOperator->getOperator(0);
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> icPelletPreconditioner(new AMP::Solver::TrilinosMLSolver(icPelletPreconditionerParams));

  boost::shared_ptr<AMP::Database> icCladPreconditioner_db = icpcSolver_db->getDatabase("cladPreconditioner"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> icCladPreconditionerParams(new AMP::Solver::SolverStrategyParameters(icCladPreconditioner_db));
  icCladPreconditionerParams->d_pOperator = columnMassOperator->getOperator(1);
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> icCladPreconditioner(new AMP::Solver::TrilinosMLSolver(icCladPreconditionerParams));

  icColumnPreconditioner->append(icPelletPreconditioner);
  icColumnPreconditioner->append(icCladPreconditioner);

  // initialize the linear solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> icLinearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(icLinearSolver_db));
  icLinearSolverParams->d_pOperator = columnMassOperator;
  icLinearSolverParams->d_comm = globalComm;
  icLinearSolverParams->d_pPreconditioner = icColumnPreconditioner;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> icLinearSolver(new AMP::Solver::PetscKrylovSolver(icLinearSolverParams));
  icLinearSolver->setZeroInitialGuess(true);
  icLinearSolver->solve(scratchVec, temperatureICPrime);
#else
  // ---------------------------------------------------------------------------------------
  // create a mass linear BVP operator
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> massICThermalModel1;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> massICThermalOp;
  boost::shared_ptr<AMP::InputDatabase> massICThermalOp_db1 =  boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MassICThermalOperator1"));
  massICThermalOp = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1, massICThermalOp_db1,massICThermalModel1));
  
    // form the rhs for the solve
  boost::shared_ptr<AMP::Database> icLinearSolver_db = input_db->getDatabase("ICLinearSolver"); 
  // ---- first initialize the preconditioner
  boost::shared_ptr<AMP::Database> icpcSolver_db = icLinearSolver_db->getDatabase("Preconditioner"); 
  boost::shared_ptr<AMP::Database> icPelletPreconditioner_db = icpcSolver_db->getDatabase("pelletPreconditioner"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> icPelletPreconditionerParams(new AMP::Solver::SolverStrategyParameters(icPelletPreconditioner_db));
  icPelletPreconditionerParams->d_pOperator = massICThermalOp;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> icPelletPreconditioner(new AMP::Solver::TrilinosMLSolver(icPelletPreconditionerParams));
  // initialize the linear solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> icLinearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(icLinearSolver_db));
  icLinearSolverParams->d_pOperator = massICThermalOp;
  icLinearSolverParams->d_comm = globalComm;
  icLinearSolverParams->d_pPreconditioner = icPelletPreconditioner;
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> icLinearSolver(new AMP::Solver::PetscKrylovSolver(icLinearSolverParams));
  icLinearSolver->setZeroInitialGuess(true);
  icLinearSolver->solve(scratchVec1, temperatureICPrime1);

#endif
  
  //-------------------------------------
  // create a  time operator for use in the preconditioner
  boost::shared_ptr<AMP::InputDatabase> timeOperator_db(new AMP::InputDatabase("TimeOperatorDatabase"));
  timeOperator_db->putDouble("CurrentDt", 0.01);
  timeOperator_db->putString("name", "TimeOperator");
  timeOperator_db->putBool("bLinearMassOperator", true);
  timeOperator_db->putBool("bLinearRhsOperator", false);
  timeOperator_db->putDouble("ScalingFactor", 1.0/0.01);

  boost::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> timeOperatorParameters(new AMP::TimeIntegrator::TimeOperatorParameters(timeOperator_db));
  timeOperatorParameters->d_pRhsOperator = linearCoupledOperator;
  timeOperatorParameters->d_pMassOperator = columnMassOperator;
  boost::shared_ptr<AMP::TimeIntegrator::TimeIntegrator::ColumnTimeOperator> columnLinearTimeOperator( new AMP::TimeIntegrator::TimeIntegrator::ColumnTimeOperator(timeOperatorParameters));
  //-------------------------------------------------------------------------//
  // initialize the column preconditioner which is a diagonal block preconditioner
  // get the ida database
  AMP_INSIST(input_db->keyExists("IDATimeIntegrator"), "Key ''IDATimeIntegrator'' is missing!");
  boost::shared_ptr<AMP::Database> ida_db = input_db->getDatabase("IDATimeIntegrator");
  boost::shared_ptr<AMP::Database> columnPreconditioner_db = ida_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(new AMP::Solver::SolverStrategyParameters(columnPreconditioner_db));
  columnPreconditionerParams->d_pOperator = columnLinearTimeOperator;
  boost::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

  boost::shared_ptr<AMP::Database> pelletPreconditioner_db = columnPreconditioner_db->getDatabase("pelletPreconditioner"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> pelletPreconditionerParams(new AMP::Solver::SolverStrategyParameters(pelletPreconditioner_db));
  pelletPreconditionerParams->d_pOperator = columnLinearTimeOperator->getOperator(0);
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> linearPelletPreconditioner(new AMP::Solver::TrilinosMLSolver(pelletPreconditionerParams));

  boost::shared_ptr<AMP::Database> cladPreconditioner_db = columnPreconditioner_db->getDatabase("cladPreconditioner"); 
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> cladPreconditionerParams(new AMP::Solver::SolverStrategyParameters(cladPreconditioner_db));
  cladPreconditionerParams->d_pOperator = columnLinearTimeOperator->getOperator(1);
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> linearCladPreconditioner(new AMP::Solver::TrilinosMLSolver(cladPreconditionerParams));

  columnPreconditioner->append(linearPelletPreconditioner);
  columnPreconditioner->append(linearCladPreconditioner);

  // ---------------------------------------------------------------------------------------
  // create the IDA time integrator
  boost::shared_ptr<AMP::TimeIntegrator::IDATimeIntegratorParameters> time_Params( new AMP::TimeIntegrator::IDATimeIntegratorParameters(ida_db));
  
  if( (time_Params.get()) == NULL ) {
    ut.failure("Testing IDATimeIntegratorParameters' Constructor");
  } else {
    ut.passes("Testing IDATimeIntegratorParameters' Constructor");
  }
    
  time_Params->d_pMassOperator = columnMassOperator;
  time_Params->d_operator = nonlinearCoupledOperator;
  time_Params->d_pPreconditioner = columnPreconditioner;
  
  time_Params->d_ic_vector = temperature;    
  time_Params->d_ic_vector_prime = temperatureICPrime;
  
  time_Params->d_pSourceTerm = RightHandSideVec;
  time_Params->d_object_name = "IDATimeIntegratorParameters";
    
  cout << "Before IDATimeIntegrator" << endl;    
  boost::shared_ptr<AMP::TimeIntegrator::IDATimeIntegrator> pIDATimeIntegrator(new AMP::TimeIntegrator::IDATimeIntegrator(time_Params));
  
  if(pIDATimeIntegrator.get() == NULL) {
    ut.failure("Testing IDATimeIntegrator's constructor");
  } else {
    ut.passes("Tested IDATimeIntegrator's constructor");
  }

  
  // ---------------------------------------------------------------------------------------
  // step in time
  int retval=0;
  double current_time=0;
  double maxT=0;
  double maxO=0;
  double abs_error=0.0;
  double minT=0;
  double minO=0;
  double rel_error=0.0;
  double exact_sol=0.0;
  int j=1;

  while(pIDATimeIntegrator->getCurrentTime() < pIDATimeIntegrator->getFinalTime())
    {
      retval = pIDATimeIntegrator->advanceSolution(pIDATimeIntegrator->getCurrentDt(), 0);
      //pIDATimeIntegrator->updateSolution();
      current_time = pIDATimeIntegrator->getCurrentTime();
      
      cout << j++ << "-th timestep" << endl;
      if(retval == 0) {
    ut.passes("Testing IDATimeIntegrator's advanceSolution. PASS!!");
      } else {
    ut.failure("Tested IDATimeIntegrator's advanceSolution. FAIL!!");
      }
      
      boost::shared_ptr<AMP::LinearAlgebra::Vector> currentSolution = pIDATimeIntegrator->getCurrentSolution();
      
      cout << "current_time = " << current_time << endl;

    }

#ifdef USE_EXT_SILO
  manager->writeFile<AMP::SiloIO> ( exeName , 0 );
#endif

  //-------------------------------------
#if 0
    if(testPassed)
  {
    ut.passes("IDA integration of thermal diffusion over clad and pellet with simple thermal contact map ");
  }
  else
  {
    ITFAILS;
  }
#endif
  input_db.reset();

  ut.passes(exeName);

//  AMP::AMPManager::shutdown();

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

  AMP::AMP_MPI::initialize();  //  PetscInitialize ( &argc , &argv , PETSC_NULL , PETSC_NULL );

  node  = AMP::AMP_MPI::getRank();
  nodes = AMP::AMP_MPI::getNodes();

    try {
    thermalContactTest(ut, argv[0]);
    gpass += ut.numPasses;
    gfail += ut.numFails;
    ut.reset();

  }
  catch (std::exception &err)
  {
    std::cout << "ERROR: " 
      << err.what()
      << std::endl;
    ut.numFails++;
  }
  catch( ... )
  {
    std::cout << "ERROR: " 
      << "An unknown exception was thrown."
      << std::endl;
    ut.numFails++;
  }

  //  PetscFinalize ();
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   



