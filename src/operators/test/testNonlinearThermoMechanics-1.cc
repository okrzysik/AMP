
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include <iostream>
#include <string>

#include "discretization/simpleDOF_Manager.h"
#include "discretization/MultiDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "operators/OperatorBuilder.h"
#include "operators/ColumnOperator.h"
#include "operators/LinearOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"

#include "operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"

#include "applyTests.h"


void thermoMechanicsTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST( input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!" );
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);

  //-----------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear mechanics
  AMP_INSIST( input_db->keyExists("testNonlinearMechanicsOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "testNonlinearMechanicsOperator", input_db, mechanicsMaterialModel));

  boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> nonlinearMechanicsVolumeOperator = boost::dynamic_pointer_cast<
    AMP::Operator::MechanicsNonlinearFEOperator>(nonlinearMechanicsOperator->getVolumeOperator());

  //---------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear thermal diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearThermalOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "testNonlinearThermalOperator", input_db, thermalTransportModel));

  //--------------------------------------------------------------------------//
  // create a column operator object for nonlinear thermomechanics
  boost::shared_ptr<AMP::Operator::OperatorParameters> params;
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearThermoMechanicsOperator(new AMP::Operator::ColumnOperator(params));
  nonlinearThermoMechanicsOperator->append(nonlinearMechanicsOperator);
  nonlinearThermoMechanicsOperator->append(nonlinearThermalOperator);

  AMP::Discretization::DOFManager::shared_ptr scalarDofMap = AMP::Discretization::simpleDOFManager::create(
      meshAdapter, AMP::Mesh::Vertex, 1, 1, true); 

  AMP::Discretization::DOFManager::shared_ptr vectorDofMap = AMP::Discretization::simpleDOFManager::create(
      meshAdapter, AMP::Mesh::Vertex, 1, 3, true); 

std::vector<AMP::Discretization::DOFManager::shared_ptr> dofMapList;
dofMapList.push_back(vectorDofMap);
dofMapList.push_back(scalarDofMap);

  AMP::Discretization::DOFManager::shared_ptr multiDofMap(
      new AMP::Discretization::multiDOFManager(meshAdapter->getComm(), dofMapList));

  // initialize the output multi-variable
  AMP::LinearAlgebra::Variable::shared_ptr multiVar = nonlinearThermoMechanicsOperator->getOutputVariable();

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector(multiDofMap, multiVar, true);
  AMP::LinearAlgebra::Vector::shared_ptr solVec = resVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = resVec->cloneVector();

  //-----------------------------------------------------------------------//
  // set up the shift and scale parameters
  double shift[2];
  double scale[2];
  shift[0] = 0.;
  shift[1] = 0.;
  scale[0] = 1.;
  scale[1] = 1.;
  std::vector<double> range(2);
  boost::shared_ptr<AMP::Operator::DiffusionTransportModel> transportModel = boost::dynamic_pointer_cast<
    AMP::Operator::DiffusionTransportModel>(thermalTransportModel);
  AMP::Materials::Material::shared_ptr matTh = transportModel->getMaterial();
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermOperator = boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearThermalOperator->getVolumeOperator());
  if ( thermOperator->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE) {
    std::string property="ThermalConductivity";
    if( (matTh->property(property))->is_argument("temperature") ) {
      range = (matTh->property(property))->get_arg_range("temperature");  // Compile error
      scale[1] = range[1]-range[0];
      shift[1] = range[0]+0.001*scale[1];
      scale[1] *= 0.999;
    }
  }

  //----------------------------------------------------------------------------//
  AMP::LinearAlgebra::Variable::shared_ptr thermVar = nonlinearThermalOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr referenceTemperatureVec = 
    AMP::LinearAlgebra::createVector(scalarDofMap, thermVar, true);
  referenceTemperatureVec->setToScalar(300.0);
  nonlinearMechanicsVolumeOperator->setReferenceTemperature(referenceTemperatureVec);

  //---------------------------------------------------------------------------//
  // now construct the linear BVP operator for mechanics
  AMP_INSIST( input_db->keyExists("testLinearMechanicsOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "testLinearMechanicsOperator", input_db, mechanicsMaterialModel));

  //--------------------------------------------------------------------------//
  // now construct the linear BVP operator for thermal
  AMP_INSIST( input_db->keyExists("testLinearThermalOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "testLinearThermalOperator", input_db, thermalTransportModel));

  //-------------------------------------------------------------------------//
  // create a column operator object for linear thermomechanics
  boost::shared_ptr<AMP::Operator::ColumnOperator> linearThermoMechanicsOperator(new AMP::Operator::ColumnOperator(params));
  linearThermoMechanicsOperator->append(linearMechanicsOperator);
  linearThermoMechanicsOperator->append(linearThermalOperator);

  ut->passes(exeName +  " : create");

  // test apply
  std::string  msgPrefix=exeName + " : apply";
  boost::shared_ptr<AMP::Operator::Operator> testOperator = nonlinearThermoMechanicsOperator;
  applyTests(ut, msgPrefix, testOperator, rhsVec, solVec, resVec, shift, scale, 2);

  ut->passes(msgPrefix);

  boost::shared_ptr<AMP::Operator::OperatorParameters> resetParams = 
    nonlinearThermoMechanicsOperator->getJacobianParameters(solVec);

  ut->passes(exeName + " : getJacobianParameters");

  linearThermoMechanicsOperator->reset(resetParams);

  std::cout<<"Tested the reset function "<<std::endl;

  ut->passes(exeName + " : Linear::reset");

}

int main(int argc, char *argv[])
{
  AMP::AMPManagerProperties startup_properties;
  startup_properties.use_MPI_Abort = false;
  AMP::AMPManager::startup(argc,argv,startup_properties);
  AMP::UnitTest ut;

  std::vector<std::string> exeNames;
  exeNames.push_back("nonlinearBVP-Mechanics-ThermalStrain-Thermal-UO2MSRZC09-1");

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


