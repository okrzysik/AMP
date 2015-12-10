#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <string>
#include "materials/Material.h"
#include "utils/shared_ptr.h"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"


#include "vectors/Variable.h"
#include "vectors/SimpleVector.h"
#include "vectors/Vector.h"

#include "utils/Writer.h"
#include "ampmesh/MeshVariable.h"

#include "operators/libmesh/MassLinearElement.h"
#include "operators/libmesh/MassLinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/MechanicsLinearFEOperator.h"
#include "operators/MechanicsNonlinearFEOperator.h"
#include "operators/MechanicsLinearElement.h"
#include "operators/MechanicsMaterialModel.h"
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
#include "operators/libmesh/VolumeIntegralOperator.h"
#include "operators/NeutronicsRhs.h"
#include "operators/CoupledOperator.h"
#include "operators/CoupledOperatorParameters.h"

#include "solvers/PetscKrylovSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"
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

  AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::PIO::logAllNodes(log_file);

  AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( mgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr pelletMeshAdapter = manager->getMesh ( "pellet" );

  AMP::LinearAlgebra::Vector::shared_ptr nullVector;

  AMP::LinearAlgebra::Variable::shared_ptr thermalVariable ( new AMP::NodalScalarVariable ( "Temperature" ) );
  AMP::LinearAlgebra::Variable::shared_ptr mechanicsVariable ( new AMP::Nodal3VectorVariable ( "Mechanics" ) );

  AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputVariable(new AMP::LinearAlgebra::MultiVariable("inputVariable"));
  inputVariable->add(thermalVariable);
  inputVariable->add(mechanicsVariable);  
  
  AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> outputVariable(new AMP::LinearAlgebra::MultiVariable("outputVariable"));
  outputVariable->add(thermalVariable);  
  outputVariable->add(mechanicsVariable);

  // first create multi-vectors for solution, time derivative, rhs, residual, and scratch 
  AMP::LinearAlgebra::Vector::shared_ptr solutionVector           = manager->createVector( inputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr solutionPrime            = manager->createVector( inputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVector                  = manager->createVector ( outputVariable );    
  AMP::LinearAlgebra::Vector::shared_ptr residualVector          = manager->createVector ( outputVariable );    
  AMP::LinearAlgebra::Vector::shared_ptr scratchVector           = manager->createVector ( inputVariable );    

  // subset the solution multivectors for ease when dealing with the individual physics
  AMP::LinearAlgebra::Vector::shared_ptr thermalSolution          = solutionVector->subsetVectorForVariable(thermalVariable);
  AMP::LinearAlgebra::Vector::shared_ptr mechanicsSolution = solutionVector->subsetVectorForVariable(mechanicsVariable);

  // subset the solution derivative multivectors for ease when dealing with the individual physics
  // note that there is no time derivative for the mechanics below
  AMP::LinearAlgebra::Vector::shared_ptr thermalPrime = solutionPrime->subsetVectorForVariable(thermalVariable);

  // subset the rhs multivectors for ease when dealing with the individual physics
  AMP::LinearAlgebra::Vector::shared_ptr thermalRhs          = rhsVector->subsetVectorForVariable(thermalVariable);
  AMP::LinearAlgebra::Vector::shared_ptr mechanicsRhs = rhsVector->subsetVectorForVariable(mechanicsVariable);
  AMP::LinearAlgebra::Vector::shared_ptr thermalScratch    = scratchVector->subsetVectorForVariable(thermalVariable); 
  
  double thermalInitialGuess = input_db->getDoubleWithDefault("ThermalInitialGuess",400);
  double referenceTemperature = input_db->getDoubleWithDefault("referenceTemperature",400);
  thermalSolution->setToScalar ( thermalInitialGuess );

  mechanicsSolution->zero();
  solutionPrime->zero();
  rhsVector->zero();
  residualVector->zero();
  scratchVector->zero();
  
  manager->registerVectorAsData ( thermalSolution , "Thermal" );
  manager->registerVectorAsData ( mechanicsSolution, "Mechanics" );
  manager->registerVectorAsData ( residualVector, "Residual" );

  //-------------------------------------
  //  create the neutronics source
  AMP_INSIST(input_db->keyExists("NeutronicsOperator"), "Key ''NeutronicsOperator'' is missing!");
  AMP::shared_ptr<AMP::Database>  neutronicsOp_db = input_db->getDatabase("NeutronicsOperator");
  AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ));
  neutronicsParams->d_MeshAdapter = pelletMeshAdapter;
  AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(new AMP::Operator::NeutronicsRhs( neutronicsParams ));

  AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar = neutronicsOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   SpecificPowerVector = manager->createVector( SpecificPowerVar );

  neutronicsOperator->apply(nullVector, nullVector, SpecificPowerVector, 1., 0.);

  //----------------------------------------------------------
  //  Integrate Nuclear Rhs over Density * Volume

  AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );

  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
  AMP::shared_ptr<AMP::Database> sourceDatabase = input_db->getDatabase("VolumeIntegralOperator");
  AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator = AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator(pelletMeshAdapter, sourceDatabase, stransportModel));

  // convert the vector of specific power to power for a given basis.
  sourceOperator->apply(nullVector, SpecificPowerVector, thermalRhs, 1., 0.);

  AMP::pout<<"L2 Norm of the thermalRhs " << thermalRhs->L2Norm()<<std::endl;

  //-----------------------------------------------
  //   create nonlinear pellet thermal operator
  //-----------------------------------------------
  
  AMP_INSIST( input_db->keyExists("pelletNonlinearThermalOperator"), "key missing!" );
  
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> pelletThermalTransportModel;
  AMP::shared_ptr<AMP::Database> pelletNonlinearThermalDatabase = input_db->getDatabase("pelletNonlinearThermalOperator");
  AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> pelletNonlinearThermalOperator = AMP::dynamic_pointer_cast<
  AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(pelletMeshAdapter, pelletNonlinearThermalDatabase, pelletThermalTransportModel));
  
  //-------------------------------------
  //   create linear pellet thermal operator
  //-------------------------------------
  
  AMP::shared_ptr<AMP::InputDatabase>  pelletLinearThermalDatabase= AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase("pelletLinearThermalOperator"));
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> pelletLinearThermalOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(pelletMeshAdapter, pelletLinearThermalDatabase, pelletThermalTransportModel));
  
  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create nonlinear pellet mechanics operator
  AMP_INSIST( input_db->keyExists("pelletNonlinearMechanicsOperator"), "key missing!" );
  
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> pelletMechanicsMaterialModel;
  AMP::shared_ptr<AMP::Database> pelletNonlinearMechanicsDatabase = input_db->getDatabase("pelletNonlinearMechanicsOperator");
  AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> pelletNonlinearMechanicsOperator = AMP::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(pelletMeshAdapter, pelletNonlinearMechanicsDatabase, pelletMechanicsMaterialModel));

  AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> pelletMechanicsVolumeOperator = AMP::dynamic_pointer_cast<
    AMP::Operator::MechanicsNonlinearFEOperator>(pelletNonlinearMechanicsOperator->getVolumeOperator());
 
  AMP::shared_ptr<AMP::LinearAlgebra::Variable> pelletBurnupVar =pelletMechanicsVolumeOperator->getInputVariable(AMP::Operator::Mechanics::BURNUP);
  AMP::shared_ptr<AMP::LinearAlgebra::Vector> pelletBurnupVector = pelletMeshAdapter->createVector(pelletBurnupVar);
  pelletBurnupVector->setToScalar(0.0001);
  pelletMechanicsVolumeOperator->setVector(AMP::Operator::Mechanics::BURNUP, pelletBurnupVector);

  AMP::LinearAlgebra::Vector::shared_ptr pelletReferenceTemperatureVec = pelletMeshAdapter->createVector( pelletMechanicsVolumeOperator->getInputVariable(AMP::Operator::Mechanics::TEMPERATURE) );
  pelletReferenceTemperatureVec->setToScalar(referenceTemperature);
  pelletMechanicsVolumeOperator->setReferenceTemperature(pelletReferenceTemperatureVec);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create linear pellet mechanics operator
  AMP_INSIST( input_db->keyExists("pelletLinearMechanicsOperator"), "key missing!" );
  AMP::shared_ptr<AMP::Database> pelletLinearMechanicsDatabase = input_db->getDatabase("pelletLinearMechanicsOperator");
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> pelletLinearMechanicsOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
      AMP::Operator::OperatorBuilder::createOperator(pelletMeshAdapter, pelletLinearMechanicsDatabase, pelletMechanicsMaterialModel));

  //-------------------------------------
  //Column of nonlinear coupled operators
  AMP::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
  AMP::shared_ptr<AMP::Operator::OperatorParameters> nonlinearParams (new AMP::Operator::OperatorParameters(tmp_db));
  AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearCoupledOperator(new AMP::Operator::ColumnOperator(nonlinearParams));
  nonlinearCoupledOperator->append(pelletNonlinearThermalOperator);
  nonlinearCoupledOperator->append(pelletNonlinearMechanicsOperator);

  //-------------------------------------
  //Column of linear operators
  AMP::shared_ptr<AMP::Operator::OperatorParameters> linearParams(new AMP::Operator::OperatorParameters(tmp_db));
  AMP::shared_ptr<AMP::Operator::ColumnOperator> linearCoupledOperator(new AMP::Operator::ColumnOperator(linearParams));
  linearCoupledOperator->append(pelletLinearThermalOperator);
  linearCoupledOperator->append(pelletLinearMechanicsOperator);

  // ---------------------------------------------------------------------------------------
  // create mass linear BVP operators
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> pelletMassThermalModel;
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> pelletMassThermalOperator;
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> pelletMassMechanicsOperator;
  
  AMP::shared_ptr<AMP::InputDatabase> massThermalOp_db1 =  AMP::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("pelletMassThermalOperator"));
  pelletMassThermalOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(pelletMeshAdapter, massThermalOp_db1,pelletMassThermalModel));
  //  pelletMassMechanicsOperator.reset();
  
  // ---------------------------------------------------------------------------------------
  // create a column mass operator object for use in the nonlinear and linear problem definition
  AMP::shared_ptr<AMP::Operator::OperatorParameters> massParams;
  AMP::shared_ptr<AMP::Operator::ColumnOperator> columnMassOperator(new AMP::Operator::ColumnOperator(massParams));
  columnMassOperator->append(pelletMassThermalOperator);
  //  columnMassOperator->append(pelletMassMechanicsOperator);
  //--------------------------------------------------------------------
  //Applying boundary conditions to the nonlinear BVP Operators

  pelletNonlinearThermalOperator->modifyRHSvector(thermalRhs);
  pelletNonlinearThermalOperator->modifyInitialSolutionVector(thermalSolution);

  pelletNonlinearMechanicsOperator->modifyRHSvector(mechanicsRhs);
  pelletNonlinearMechanicsOperator->modifyInitialSolutionVector(mechanicsSolution);

  // call apply, it serves two purposes, 1) called before the getJacobian for mechanics(not sure if needed)
  //                                                       2)  form the rhs for the solve for the time derivative of temperature
  nonlinearCoupledOperator->apply(rhsVector, solutionVector, scratchVector, -1.0, 1.0);
  // unfortunately we have to apply a 'trick' here to zero out the Dirichlet boundaries
  AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletBoundaryOperator = AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(pelletNonlinearThermalOperator->getBoundaryOperator());
  //  AMP_INSIST(dirichletBoundaryOperator.get()!=NULL);
  dirichletBoundaryOperator->apply(nullVector, nullVector, scratchVector, 1.0, 0.0);
  
  pelletLinearMechanicsOperator->reset(pelletNonlinearMechanicsOperator->getParameters("Jacobian", mechanicsSolution));
  //-------------------------------------------------------------------------------------------
  // solve for the time derivative for thermal at the initial time

  // create a mass linear BVP operator
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> pelletMassICThermalModel;
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> pelletMassICThermalOperator;
  AMP::shared_ptr<AMP::InputDatabase> pelletMassICThermalOp_db =  AMP::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("pelletMassICThermalOperator"));
  pelletMassICThermalOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(pelletMeshAdapter, pelletMassICThermalOp_db,pelletMassICThermalModel));
  
    // form the rhs for the solve
  AMP::shared_ptr<AMP::Database> icLinearSolver_db = input_db->getDatabase("ICLinearSolver"); 
  // ---- first initialize the preconditioner
  AMP::shared_ptr<AMP::Database> icpcSolver_db = icLinearSolver_db->getDatabase("Preconditioner"); 
  AMP::shared_ptr<AMP::Database> icPelletPreconditioner_db = icpcSolver_db->getDatabase("pelletPreconditioner"); 
  AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> icPelletPreconditionerParams(new AMP::Solver::SolverStrategyParameters(icPelletPreconditioner_db));
  icPelletPreconditionerParams->d_pOperator = pelletMassICThermalOperator;
  AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> icPelletPreconditioner(new AMP::Solver::TrilinosMLSolver(icPelletPreconditionerParams));
  // initialize the linear solver
  AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> icLinearSolverParams(new
      AMP::Solver::PetscKrylovSolverParameters(icLinearSolver_db));
  icLinearSolverParams->d_pOperator = pelletMassICThermalOperator;
  icLinearSolverParams->d_comm = globalComm;
  icLinearSolverParams->d_pPreconditioner = icPelletPreconditioner;
  AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> icLinearSolver(new AMP::Solver::PetscKrylovSolver(icLinearSolverParams));
  icLinearSolver->setZeroInitialGuess(true);
  icLinearSolver->solve(thermalScratch, thermalPrime);

  //-------------------------------------
  // create a  time operator for use in the preconditioner
  AMP_INSIST(input_db->keyExists("IDATimeIntegrator"), "Key ''IDATimeIntegrator'' is missing!");
  AMP::shared_ptr<AMP::Database> ida_db = input_db->getDatabase("IDATimeIntegrator");
  AMP::shared_ptr<AMP::InputDatabase> timeOperator_db(new AMP::InputDatabase("TimeOperatorDatabase"));
  timeOperator_db->putDouble("CurrentDt", 0.01);
  timeOperator_db->putString("name", "TimeOperator");
  timeOperator_db->putBool("bLinearMassOperator", true);
  timeOperator_db->putBool("bLinearRhsOperator", false);
  timeOperator_db->putDouble("ScalingFactor", 1.0/0.01);
  timeOperator_db->putInteger("algebraicComponent", ida_db->getIntegerWithDefault("algebraicComponent", -1));

  AMP::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> timeOperatorParameters(new AMP::TimeIntegrator::TimeOperatorParameters(timeOperator_db));
  timeOperatorParameters->d_pRhsOperator = linearCoupledOperator;
  timeOperatorParameters->d_pMassOperator = columnMassOperator;
  AMP::shared_ptr<AMP::TimeIntegrator::TimeIntegrator::ColumnTimeOperator> columnLinearTimeOperator( new AMP::TimeIntegrator::TimeIntegrator::ColumnTimeOperator(timeOperatorParameters));
  //-------------------------------------------------------------------------//
  // initialize the column preconditioner which is a diagonal block preconditioner
  // get the ida database
  AMP::shared_ptr<AMP::Database> columnPreconditioner_db = ida_db->getDatabase("Preconditioner");
  AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(new AMP::Solver::SolverStrategyParameters(columnPreconditioner_db));
  columnPreconditionerParams->d_pOperator = columnLinearTimeOperator;
  AMP::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

  AMP::shared_ptr<AMP::Database> pelletThermalPreconditioner_db = columnPreconditioner_db->getDatabase("pelletThermalPreconditioner"); 
  AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> pelletThermalPreconditionerParams(new AMP::Solver::SolverStrategyParameters(pelletThermalPreconditioner_db));
  pelletThermalPreconditionerParams->d_pOperator = columnLinearTimeOperator->getOperator(0);
  AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> pelletLinearThermalPreconditioner(new AMP::Solver::TrilinosMLSolver(pelletThermalPreconditionerParams));

  AMP::shared_ptr<AMP::Database> pelletMechanicsPreconditioner_db = columnPreconditioner_db->getDatabase("pelletMechanicsPreconditioner"); 
  AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> pelletMechanicsPreconditionerParams(new AMP::Solver::SolverStrategyParameters(pelletMechanicsPreconditioner_db));
  pelletMechanicsPreconditionerParams->d_pOperator = columnLinearTimeOperator->getOperator(1);
  AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> pelletLinearMechanicsPreconditioner(new AMP::Solver::TrilinosMLSolver(pelletMechanicsPreconditionerParams));

  columnPreconditioner->append(pelletLinearThermalPreconditioner);
  columnPreconditioner->append(pelletLinearMechanicsPreconditioner);

  // ---------------------------------------------------------------------------------------
  // create the IDA time integrator
  AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegratorParameters> time_Params( new AMP::TimeIntegrator::IDATimeIntegratorParameters(ida_db));
  
  if( (time_Params.get()) == NULL ) {
    ut.failure("Testing IDATimeIntegratorParameters' Constructor");
  } else {
    ut.passes("Testing IDATimeIntegratorParameters' Constructor");
  }
    
  time_Params->d_pMassOperator = columnMassOperator;
  time_Params->d_operator = nonlinearCoupledOperator;
  time_Params->d_pPreconditioner = columnPreconditioner;
  
  time_Params->d_ic_vector = solutionVector;    
  time_Params->d_ic_vector_prime = solutionPrime;
  
  time_Params->d_pSourceTerm = rhsVector;
  time_Params->d_object_name = "IDATimeIntegratorParameters";
  time_Params->d_pAlgebraicVariable = pelletNonlinearMechanicsOperator->getOutputVariable();
  
  cout << "Before IDATimeIntegrator" << endl;    
  AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegrator> pIDATimeIntegrator(new AMP::TimeIntegrator::IDATimeIntegrator(time_Params));
  
  if(pIDATimeIntegrator.get() == NULL) {
    ut.failure("Testing IDATimeIntegrator's constructor");
  } else {
    ut.passes("Tested IDATimeIntegrator's constructor");
  }

  AMP::shared_ptr<AMP::Operator::Operator> pIDAOperator = pIDATimeIntegrator->getOperator();

  pIDAOperator->apply(nullVector, solutionVector, scratchVector);
  AMP::pout << "Residual at initial time is " << scratchVector->L2Norm() << std::endl;
  
  // ---------------------------------------------------------------------------------------
  // step in time
  int retval=0;
  double current_time=0;
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
      
      AMP::shared_ptr<AMP::LinearAlgebra::Vector> currentSolution = pIDATimeIntegrator->getCurrentSolution();
      
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

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   



