#include <string>
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "materials/Material.h"
#include "utils/shared_ptr.h"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "operators/NeutronicsRhs.h"
#include "vectors/Variable.h"

#include "utils/Writer.h"
#include "vectors/Vector.h"
#include "operators/libmesh/VolumeIntegralOperator.h"
#include "operators/NeutronicsRhs.h"
#include "operators/NonlinearBVPOperator.h"

/* libMesh files */
#include "libmesh.h"
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

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/libmesh/MassLinearFEOperator.h"

#include "operators/IsotropicElasticModel.h"
#include "operators/MechanicsLinearElement.h"
#include "operators/MechanicsLinearFEOperator.h"
#include "operators/DirichletMatrixCorrection.h"
#include "operators/DirichletVectorCorrection.h"
#include "operators/LinearBVPOperator.h"
#include "operators/ColumnOperator.h"
#include "solvers/ColumnSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"
#include "time_integrators/ColumnTimeOperator.h"
#include "time_integrators/sundials/IDATimeIntegrator.h"
#include "time_integrators/sundials/IDATimeOperator.h"

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

#define __PI__ 3.14159265

//#define __INIT_FN__(x, y, z, t) ( exp(- __PI__ * __PI__ * t) * sin(__PI__ * x) * sin(__PI__ * y) * sin(__PI__ * z) )
//#define __INIT_FN__(x,y,z,t) ( exp(-0.015 *  __PI__ * __PI__ * t) * cos(0.1 * __PI__ * x) * cos(0.1 * __PI__ * y) * cos(0.05 * __PI__ * z) )
#define __INIT_FN__(x,y,z,t) ( 750.0+ 10000.0*(0.5+ x) * (0.5 -x) *(0.5+ y) * (0.5 -y) *(0.5+ z) * (0.5 -z) )

void IDATimeIntegratorTest(AMP::UnitTest *ut )
{
  std::string input_file = "input_testIDA-NonlinearColumnWithTwoIndependentOperators-1";
  std::string log_file = "output_testIDA-NonlinearColumnWithTwoIndependentOperators-1";

  AMP::PIO::logOnlyNodeZero(log_file);

  AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "ida" );
  
  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create two nonlinear thermal BVP operators
  AMP_INSIST( input_db->keyExists("NonlinearThermalOperator1"), "key missing!" );
  AMP_INSIST( input_db->keyExists("NonlinearThermalOperator2"), "key missing!" );

  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel1;
  AMP::shared_ptr<AMP::Database> nonlinearDatabase1 = input_db->getDatabase("NonlinearThermalOperator1");
  AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator1 = AMP::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, nonlinearDatabase1, elementModel1));
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel2;
  AMP::shared_ptr<AMP::Database> nonlinearDatabase2 = input_db->getDatabase("NonlinearThermalOperator2");
  AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator2 = AMP::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, nonlinearDatabase2, elementModel2));
  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a column rhs operator object with the nonlinear thermal in it for use in the nonlinear problem definition
  AMP::shared_ptr<AMP::Operator::OperatorParameters> params;
  AMP::shared_ptr<AMP::Operator::ColumnOperator> columnNonlinearRhsOperator(new AMP::Operator::ColumnOperator(params));
  columnNonlinearRhsOperator->append(nonlinearThermalOperator1);
  columnNonlinearRhsOperator->append(nonlinearThermalOperator2);
  // ---------------------------------------------------------------------------------------
  // create linear BVP operators
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator1;
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator2;
  AMP::shared_ptr<AMP::InputDatabase> linearOp1_db =  AMP::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("LinearThermalOperator1"));
  AMP::shared_ptr<AMP::InputDatabase> linearOp2_db =  AMP::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("LinearThermalOperator2"));
  linearThermalOperator1 = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, linearOp1_db, elementModel1));
  linearThermalOperator2 = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, linearOp2_db, elementModel2));
  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a column rhs operator object with the linear thermal in it for use in the linear problem definition
  AMP::shared_ptr<AMP::Operator::OperatorParameters> rhsparams;
  AMP::shared_ptr<AMP::Operator::ColumnOperator> columnLinearRhsOperator(new AMP::Operator::ColumnOperator(rhsparams));
  columnLinearRhsOperator->append(linearThermalOperator1);
  columnLinearRhsOperator->append(linearThermalOperator2);
  
  // ---------------------------------------------------------------------------------------
  // create a mass linear BVP operator
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> massElementModel1;
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> massElementModel2;
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> massOperator1;
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> massOperator2;
  AMP::shared_ptr<AMP::InputDatabase> massOp1_db =  AMP::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MassOperator1"));
  AMP::shared_ptr<AMP::InputDatabase> massOp2_db =  AMP::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MassOperator2"));
  massOperator1 = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, massOp1_db,massElementModel1));
  massOperator2 = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, massOp2_db,massElementModel2));
  
  // ---------------------------------------------------------------------------------------
  // create a column mass operator object for use in the nonlinear and linear problem definition
  AMP::shared_ptr<AMP::Operator::OperatorParameters> massParams;
  AMP::shared_ptr<AMP::Operator::ColumnOperator> columnMassOperator(new AMP::Operator::ColumnOperator(massParams));
  columnMassOperator->append(massOperator1);
  columnMassOperator->append(massOperator2);
  
  // ---------------------------------------------------------------------------------------
  // create a  time operator for use in the preconditioner
  AMP::shared_ptr<AMP::InputDatabase> timeOperator_db(new AMP::InputDatabase("TimeOperatorDatabase"));
  timeOperator_db->putDouble("CurrentDt", 0.01);
  timeOperator_db->putString("name", "TimeOperator");
  timeOperator_db->putBool("bLinearMassOperator", true);
  timeOperator_db->putBool("bLinearRhsOperator", false);
  timeOperator_db->putDouble("ScalingFactor", 1.0/0.01);

  AMP::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> timeOperatorParameters(new AMP::TimeIntegrator::TimeOperatorParameters(timeOperator_db));
  timeOperatorParameters->d_pRhsOperator = columnLinearRhsOperator;
  timeOperatorParameters->d_pMassOperator = columnMassOperator;
  AMP::shared_ptr<AMP::TimeIntegrator::TimeIntegrator::ColumnTimeOperator> columnLinearTimeOperator( new AMP::TimeIntegrator::TimeIntegrator::ColumnTimeOperator(timeOperatorParameters));
  
  // ---------------------------------------------------------------------------------------
  // create vectors for initial conditions (IC) and time derivative at IC
  AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator1 = AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearThermalOperator1->getVolumeOperator());
  AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator2 = AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearThermalOperator2->getVolumeOperator());

  // note that the input variable for the time integrator and time operator will be a multivariable
  AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputVar(new AMP::LinearAlgebra::MultiVariable("inputVariable"));
  inputVar->add(thermalVolumeOperator1->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
  inputVar->add(thermalVolumeOperator2->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
  
  AMP::LinearAlgebra::Variable::shared_ptr outputVar = columnNonlinearRhsOperator->getOutputVariable();
  
  AMP::LinearAlgebra::Vector::shared_ptr  initialCondition = meshAdapter->createVector(inputVar);
  AMP::LinearAlgebra::Vector::shared_ptr  initialConditionPrime = meshAdapter->createVector(inputVar);
  AMP::LinearAlgebra::Vector::shared_ptr   f = meshAdapter->createVector(outputVar);

  // ---------------------------------------------------------------------------------------
  //  create neutronics source
  AMP_INSIST(input_db->keyExists("NeutronicsOperator"), "Key ''NeutronicsOperator'' is missing!");
  AMP::shared_ptr<AMP::Database>  neutronicsOp_db = input_db->getDatabase("NeutronicsOperator");
  AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ));
  neutronicsParams->d_MeshAdapter = meshAdapter;
  AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(new AMP::Operator::NeutronicsRhs( neutronicsParams ));
  
  AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar = neutronicsOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr   SpecificPowerVec = meshAdapter->createVector( SpecificPowerVar );

  // create the following shared pointers for ease of use
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  neutronicsOperator->apply( nullVec, SpecificPowerVec);
  
  //  Integrate Nuclear Rhs over Density * Volume //
  
  AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );
  
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> sourceTransportModel;
  AMP::shared_ptr<AMP::Database> sourceDatabase = input_db->getDatabase("VolumeIntegralOperator");
  AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator = AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, sourceDatabase, sourceTransportModel));
  
  // Create the power (heat source) vector.
  AMP::LinearAlgebra::Variable::shared_ptr powerInWattsVar = sourceOperator->createOutputVariable("HeatSource");
  AMP::LinearAlgebra::Vector::shared_ptr   powerInWattsVec = meshAdapter->createVector( powerInWattsVar );
  powerInWattsVec->zero();
  
  // convert the vector of specific power to power for a given basis.
  sourceOperator->apply(SpecificPowerVec, powerInWattsVec);
  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // set initial conditions, initialize created vectors

  AMP::Mesh::MeshManager::Adapter::NodeIterator node = meshAdapter->beginNode();
  AMP::Mesh::MeshManager::Adapter::NodeIterator end_node = meshAdapter->endNode();
  
  AMP::Mesh::DOFMap::shared_ptr dof_map = meshAdapter->getDOFMap(thermalVolumeOperator1->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));

  AMP::shared_ptr<AMP::LinearAlgebra::Vector> thermalIC1 = initialCondition->subsetVectorForVariable(thermalVolumeOperator1->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
  AMP::shared_ptr<AMP::LinearAlgebra::Vector> thermalIC2 = initialCondition->subsetVectorForVariable(thermalVolumeOperator2->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
  int counter=0;     
  for( ; node != end_node ; ++node)
    {
      counter++;
      
      std::vector<unsigned int> bndGlobalIds;
      std::vector<unsigned int> d_dofIds;
      d_dofIds.resize(0);
      dof_map->getDOFs(*node, bndGlobalIds, d_dofIds);
            
      double px = (*node).x();
      double py = (*node).y();
      double pz = (*node).z();
      
      double val = __INIT_FN__(px, py, pz, 0);
      cout << "val = " << val << endl;
      
      cout << "counter = " << counter << "bndGlobalIds.size() = " << bndGlobalIds.size() << endl;
      for(unsigned int i = 0; i < bndGlobalIds.size(); i++)
    {
      thermalIC1->setValueByGlobalID(bndGlobalIds[i], val);
      thermalIC2->setValueByGlobalID(bndGlobalIds[i], val);
    }//end for i
    }//end for node
  
  // ** please do not set the time derivative to be non-zero!!
  // ** as this causes trouble with the boundary - BP, 07/16/2010
  initialConditionPrime->zero();

  AMP::shared_ptr<AMP::LinearAlgebra::Vector> thermalRhs1 = f->subsetVectorForVariable(thermalVolumeOperator1->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
  AMP::shared_ptr<AMP::LinearAlgebra::Vector> thermalRhs2 = f->subsetVectorForVariable(thermalVolumeOperator2->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
  // create a copy of the rhs which can be modified at each time step (maybe)
  thermalRhs1->copyVector(powerInWattsVec);
  thermalRhs2->copyVector(powerInWattsVec);
  // modify the rhs to take into account boundary conditions
  nonlinearThermalOperator1->modifyRHSvector(f);
  nonlinearThermalOperator1->modifyInitialSolutionVector(initialCondition);
  nonlinearThermalOperator2->modifyRHSvector(f);
  nonlinearThermalOperator2->modifyInitialSolutionVector(initialCondition);
  
  // ---------------------------------------------------------------------------------------
  // create a preconditioner

  // get the ida database
  AMP_INSIST(input_db->keyExists("IDATimeIntegrator"), "Key ''IDATimeIntegrator'' is missing!");
  AMP::shared_ptr<AMP::Database> ida_db = input_db->getDatabase("IDATimeIntegrator");
  // initialize the column preconditioner which is a diagonal block preconditioner
  AMP::shared_ptr<AMP::Database> columnPreconditioner_db = ida_db->getDatabase("Preconditioner");
  AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(new AMP::Solver::SolverStrategyParameters(columnPreconditioner_db));
  if(columnPreconditionerParams.get() == NULL) {
    ut.failure("Testing SolverStrategyParameters's constructor: FAIL");
  } else {
    ut.passes("Testing SolverStrategyParameters's constructor: PASS");
  }
  
  columnPreconditionerParams->d_pOperator = columnLinearTimeOperator;
  AMP::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

  AMP::shared_ptr<AMP::Database> thermalPreconditioner1_db = columnPreconditioner_db->getDatabase("thermalPreconditioner1"); 
  AMP::shared_ptr<AMP::Database> thermalPreconditioner2_db = columnPreconditioner_db->getDatabase("thermalPreconditioner2"); 
  AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams1(new AMP::Solver::SolverStrategyParameters(thermalPreconditioner1_db));
  AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams2(new AMP::Solver::SolverStrategyParameters(thermalPreconditioner2_db));
  thermalPreconditionerParams1->d_pOperator = columnLinearTimeOperator->getOperator(0);
  thermalPreconditionerParams2->d_pOperator = columnLinearTimeOperator->getOperator(1);
  AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalPreconditioner1(new AMP::Solver::TrilinosMLSolver(thermalPreconditionerParams1));
  AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalPreconditioner2(new AMP::Solver::TrilinosMLSolver(thermalPreconditionerParams2));

  columnPreconditioner->append(linearThermalPreconditioner1);
  columnPreconditioner->append(linearThermalPreconditioner2);
  
  if(columnPreconditioner.get() == NULL) {
    ut.failure("Testing column preconditioner's constructor: FAIL");
  } else {
    ut.passes("Testing column preconditioner's constructor: PASS");
  }

  // ---------------------------------------------------------------------------------------
  // create the IDA time integrator
  AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegratorParameters> time_Params( new AMP::TimeIntegrator::IDATimeIntegratorParameters(ida_db));
  
  if( (time_Params.get()) == NULL ) {
    ut.failure("Testing IDATimeIntegratorParameters' Constructor");
  } else {
    ut.passes("Testing IDATimeIntegratorParameters' Constructor");
  }
    
  time_Params->d_pMassOperator = columnMassOperator;
  time_Params->d_operator = columnNonlinearRhsOperator;
  time_Params->d_pPreconditioner = columnPreconditioner;
  
  time_Params->d_ic_vector = initialCondition;    
  time_Params->d_ic_vector_prime = initialConditionPrime;
  
  time_Params->d_pSourceTerm = f;
  time_Params->d_object_name = "IDATimeIntegratorParameters";
    
  cout << "Before IDATimeIntegrator" << endl;    
  AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegrator> pIDATimeIntegrator(new AMP::TimeIntegrator::IDATimeIntegrator(time_Params));
  
  if(pIDATimeIntegrator.get() == NULL) {
    ut.failure("Testing IDATimeIntegrator's constructor");
  } else {
    ut.passes("Tested IDATimeIntegrator's constructor");
  }
  
  // ---------------------------------------------------------------------------------------
  // step in time
  int retval=0;
  double current_time=0;
  double max=0;
  double min=0;
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
      
      max = pIDATimeIntegrator->getCurrentSolution()->max();
      min = pIDATimeIntegrator->getCurrentSolution()->min();
      
      cout << "current_time = " << current_time << endl;
      cout << "max val of the current solution = " << max << endl;
      cout << "min val of the current solution = " << min << endl;
    }
  
  
  AMP::AMPManager::shutdown();
  
  if (ut.numFails == 0)
    {
      ut.passes("testIDATimeIntegrator successful");
    }
}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {        
        IDATimeIntegratorTest(ut);        
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

//---------------------------------------------------------------------------//
//                        end of SundialsVectorTest.cc
//---------------------------------------------------------------------------//








