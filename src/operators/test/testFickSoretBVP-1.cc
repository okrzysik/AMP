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

#include "ampmesh/MeshVariable.h"

#include "libmesh.h"

#include "operators/diffusion/FickSoretNonlinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/NeumannVectorCorrection.h"

#include "../BVPOperatorParameters.h"
#include "../LinearBVPOperator.h"
#include "../NonlinearBVPOperator.h"
#include "../OperatorBuilder.h"

#include "applyTests.h"


void bvpTest1(AMP::UnitTest *ut, const std::string exeName)
{
  // Tests FickSoret Dirchlet BVP operator for temperature

  // Initialization
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  // Input database
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  std::string mesh_file = input_db->getString("Mesh");

  // Mesh
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = AMP::Mesh::MeshManager::Adapter::shared_ptr ( new AMP::Mesh::MeshManager::Adapter () );
  meshAdapter->readExodusIIFile ( mesh_file.c_str() );

  // Create nonlinear FickSoret BVP operator and access volume nonlinear FickSoret operator
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::Operator> nlinBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter, "testFickSoretBVPOperator", input_db, elementPhysicsModel);
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nlinBVPOp =
          boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(nlinBVPOperator);
  boost::shared_ptr<AMP::Operator::FickSoretNonlinearFEOperator> nlinOp =
         boost::dynamic_pointer_cast<AMP::Operator::FickSoretNonlinearFEOperator>(nlinBVPOp->getVolumeOperator());
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickOp =
         boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinOp->getFickOperator());
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> soretOp =
         boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinOp->getSoretOperator());
  
  // use the linear BVP operator to create a Fick linear operator with bc's
  boost::shared_ptr<AMP::Operator::Operator> linBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter, "testLinearFickBVPOperator", input_db, elementPhysicsModel);
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linBVPOp =
          boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(linBVPOperator);

  ut->passes(exeName+": creation");
  std::cout.flush();

  // set shift, scale for applyTests
  double shift[2];
  double scale[2];
  shift[0] = 0.;
  shift[1] = 0.;
  scale[0] = 1.;
  scale[1] = 1.;
  std::vector<double> trange(2), crange(2);
  boost::shared_ptr<AMP::Operator::DiffusionTransportModel> fickTransportModel = fickOp->getTransportModel();
  std::vector<double> defaults(2,0);
  AMP::Materials::Material::shared_ptr fmat = fickTransportModel->getMaterial();
  // the Soret has a principal variable of temperature
  if (soretOp->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE) {
      std::string property="ThermalDiffusionCoefficient";
      if( (fmat->property(property))->is_argument("temperature") ) {
          trange = (fmat->property(property))->get_arg_range("temperature");  // Compile error
          scale[1] = trange[1]-trange[0];
          shift[1] = trange[0]+0.001*scale[1];
          scale[1] *= 0.999;
          defaults = (fmat->property(property))->get_defaults();  // compile error
      }
  }
  // the Fick has a principal variable of temperature
  if (fickOp->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION) {
      std::string property="FickCoefficient";
      if( (fmat->property(property))->is_argument("concentration") ) {
          crange = (fmat->property(property))->get_arg_range("concentration");  // Compile error
          scale[0] = crange[1]-crange[0];
          shift[0] = crange[0]+0.001*scale[0];
          scale[0] *= 0.999;
          defaults = (fmat->property(property))->get_defaults();  // compile error
      }
  }

  // Set up input and output vectors
  AMP::LinearAlgebra::Variable::shared_ptr cVar(fickOp->getInputVariable(AMP::Operator::Diffusion::CONCENTRATION));
  AMP::LinearAlgebra::Variable::shared_ptr tVar(soretOp->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
  boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> fsInpVar(new AMP::LinearAlgebra::MultiVariable("fsInput"));
  fsInpVar->add(cVar);
  fsInpVar->add(tVar);
  boost::shared_ptr<AMP::LinearAlgebra::Variable> fsOutVar(nlinBVPOp->getOutputVariable());

  AMP::LinearAlgebra::Vector::shared_ptr solVec = meshAdapter->createVector(fsInpVar);
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = meshAdapter->createVector(fsOutVar);
  AMP::LinearAlgebra::Vector::shared_ptr resVec = meshAdapter->createVector(fsOutVar);

  // set default values of input variables
  AMP::LinearAlgebra::Vector::shared_ptr inConcVec = solVec->subsetVectorForVariable(cVar);
  AMP::LinearAlgebra::Vector::shared_ptr inTempVec = solVec->subsetVectorForVariable(tVar);
  inConcVec->setToScalar(defaults[1]);  // compile error
  inTempVec->setToScalar(defaults[0]);  // compile error
  rhsVec->setToScalar(0.);

  AMP_INSIST(nlinOp->isValidInput(solVec), "input variable not set up correctly");

  // Test apply
  std::string msgPrefix=exeName+": apply";
  {
    applyTests(ut, msgPrefix, nlinBVPOperator, rhsVec, solVec, resVec, shift, scale, 2);
  }
  std::cout.flush();

  // Test linear reset from getJacobianParameters
  for(int i = 0; i < 3; i++) {
    inConcVec->setRandomValues();
    inTempVec->setRandomValues();
    adjust(solVec, shift, scale, 2);
    rhsVec->setRandomValues();
    resVec->setRandomValues();
    boost::shared_ptr<AMP::Operator::OperatorParameters> jacparams = nlinBVPOp->getJacobianParameters(solVec);
    linBVPOp->reset(jacparams);
  }//end for i

  ut->passes(exeName+": getJacobianParameters");
  std::cout.flush();

}

int main(int argc, char *argv[])
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);

    AMP::UnitTest ut;

  const int NUMFILES=2;
  std::string files[NUMFILES] = {
        "FickSoret-BVP-TUI-1", "FickSoret-BVP-UO2MSRZC09-1"
  };

    try {
    for (int i=0; i<NUMFILES; i++) {
        bvpTest1(&ut, files[i]);    }

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


