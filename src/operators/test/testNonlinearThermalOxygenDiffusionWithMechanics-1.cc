#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshVariable.h"

#include "materials/Material.h"

#include "operators/OperatorBuilder.h"
#include "operators/ColumnOperator.h"
#include "operators/LinearOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/LinearBVPOperator.h"

#include "operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"

#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"

#include "applyTests.h"


void thermoMechanicsTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "cylinder" );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear mechanics
  AMP_INSIST( input_db->keyExists("testNonlinearMechanicsOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														    "testNonlinearMechanicsOperator",
														    input_db,
														    mechanicsMaterialModel));


  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear thermal diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearThermalOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														    "testNonlinearThermalOperator",
														    input_db,
														    thermalTransportModel));
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermOperator = boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearThermalOperator->getVolumeOperator());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear oxygen diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearOxygenOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> oxygenTransportModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearOxygenOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														    "testNonlinearOxygenOperator",
														    input_db,
														    oxygenTransportModel));
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickOperator = boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearOxygenOperator->getVolumeOperator());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a column operator object for nonlinear thermal and oxygen diffusion, and mechanics
  boost::shared_ptr<AMP::Operator::OperatorParameters> params;
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearThermalOxygenDiffusionMechanicsOperator(new AMP::Operator::ColumnOperator(params));
  nonlinearThermalOxygenDiffusionMechanicsOperator->append(nonlinearMechanicsOperator);
  nonlinearThermalOxygenDiffusionMechanicsOperator->append(nonlinearThermalOperator);
  nonlinearThermalOxygenDiffusionMechanicsOperator->append(nonlinearOxygenOperator);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the input multi-variable
  boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> volumeOperator = boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(nonlinearMechanicsOperator->getVolumeOperator());
  boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputVariable(new AMP::LinearAlgebra::MultiVariable("inputVariable"));
  inputVariable->add(volumeOperator->getInputVariable(AMP::Operator::Mechanics::DISPLACEMENT));
  inputVariable->add(volumeOperator->getInputVariable(AMP::Operator::Mechanics::TEMPERATURE));
  inputVariable->add(volumeOperator->getInputVariable(AMP::Operator::Mechanics::OXYGEN_CONCENTRATION));

  // initialize the output multi-variable
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable = nonlinearThermalOxygenDiffusionMechanicsOperator->getOutputVariable();

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr solVec = meshAdapter->createVector( inputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = meshAdapter->createVector( outputVariable );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = meshAdapter->createVector( outputVariable );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // set up the frozen variables for each operator
  // first get defaults
  double defTemp, defConc;
  boost::shared_ptr<AMP::Operator::DiffusionTransportModel> transportModelTh =
    boost::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>(thermalTransportModel);
  defConc = transportModelTh->getDefault(AMP::Operator::Diffusion::CONCENTRATION);
  boost::shared_ptr<AMP::Operator::DiffusionTransportModel> transportModelOx =
    boost::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>(oxygenTransportModel);
  defTemp = transportModelOx->getDefault(AMP::Operator::Diffusion::TEMPERATURE);
  // next get vectors
  AMP::LinearAlgebra::Vector::shared_ptr tempVec = solVec->subsetVectorForVariable(inputVariable->getVariable(1));
  AMP::LinearAlgebra::Vector::shared_ptr concVec = solVec->subsetVectorForVariable(inputVariable->getVariable(2));
  tempVec->setToScalar(defTemp);
  concVec->setToScalar(defConc);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // set up the shift and scale parameters
  double shift[2];
  double scale[2];
  shift[0] = 0.;
  shift[1] = 0.;
  scale[0] = 1.;
  scale[1] = 1.;
  std::vector<double> range(2);
  AMP::Materials::Material::shared_ptr matTh = transportModelTh->getMaterial();
  AMP::Materials::Material::shared_ptr matOx = transportModelOx->getMaterial();
  if ( thermOperator->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE) {
      std::string property="ThermalConductivity";
      if( (matTh->property(property))->is_argument("temperature") ) {
          range = (matTh->property(property))->get_arg_range("temperature");  // Compile error
          scale[1] = range[1]-range[0];
          shift[1] = range[0]+0.001*scale[1];
          scale[1] *= 0.999;
      }
  }
  // the Fick has a principal variable of oxygen
  if ( fickOperator->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION) {
      std::string property="FickCoefficient";
      if( (matOx->property(property))->is_argument("concentration") ) {
          range = (matOx->property(property))->get_arg_range("concentration");  // Compile error
          scale[0] = range[1]-range[0];
          shift[0] = range[0]+0.001*scale[0];
          scale[0] *= 0.999;
      }
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // IMPORTANT:: call init before proceeding any further on the nonlinear mechanics operator
  AMP::LinearAlgebra::Vector::shared_ptr referenceTemperatureVec = meshAdapter->createVector( volumeOperator->getInputVariable(AMP::Operator::Mechanics::TEMPERATURE) );
  referenceTemperatureVec->setToScalar(300.0);
  volumeOperator->setReferenceTemperature(referenceTemperatureVec);
  volumeOperator->init();
  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for mechanics
  AMP_INSIST( input_db->keyExists("testLinearMechanicsOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														 "testLinearMechanicsOperator",
														 input_db,
														 mechanicsMaterialModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for thermal
  AMP_INSIST( input_db->keyExists("testLinearThermalOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														 "testLinearThermalOperator",
														 input_db,
														 thermalTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for oxygen
  AMP_INSIST( input_db->keyExists("testLinearOxygenOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearOxygenOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
														 "testLinearOxygenOperator",
														 input_db,
														 oxygenTransportModel));
  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a column operator object for linear thermomechanics
  boost::shared_ptr<AMP::Operator::ColumnOperator> linearThermalOxygenDiffusionMechanicsOperator(new AMP::Operator::ColumnOperator(params));
  linearThermalOxygenDiffusionMechanicsOperator->append(linearMechanicsOperator);
  linearThermalOxygenDiffusionMechanicsOperator->append(linearThermalOperator);
  linearThermalOxygenDiffusionMechanicsOperator->append(linearOxygenOperator);

  ut->passes(exeName +  " : create");

  // test apply
  std::string  msgPrefix=exeName + " : apply";
  boost::shared_ptr<AMP::Operator::Operator> testOperator = boost::dynamic_pointer_cast<AMP::Operator::Operator>(nonlinearThermalOxygenDiffusionMechanicsOperator);
  applyTests(ut, msgPrefix, testOperator, rhsVec, solVec, resVec, shift, scale, 3);

  ut->passes(msgPrefix);

  boost::shared_ptr<AMP::Operator::OperatorParameters> resetParams = nonlinearThermalOxygenDiffusionMechanicsOperator->getJacobianParameters(solVec);

  ut->passes(exeName + " : getJacobianParameters");

  linearThermalOxygenDiffusionMechanicsOperator->reset(resetParams);

  ut->passes(exeName + " : Linear::reset");

}

int main(int argc, char *argv[])
{
  AMP::AMPManagerProperties startup_properties;
  startup_properties.use_MPI_Abort = false;
  AMP::AMPManager::startup(argc,argv,startup_properties);
  AMP::UnitTest ut;

  std::vector<std::string> exeNames;
  exeNames.push_back("nonlinearBVP-Mechanics-ThermalStrain-Thermal-Oxygen-UO2MSRZC09-1");
  //  exeNames.push_back("testNonlinearMechanics-1-reduced");

  for(unsigned int i = 0; i < exeNames.size(); i++) {
    try {
      thermoMechanicsTest(&ut, exeNames[i]);
    } catch (std::exception &err) {
      std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
      ut.failure("ERROR: While testing");
    } catch( ... ) {
      std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
      ut.failure("ERROR: While testing");
    }
  }

  ut.report();

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;
}   


