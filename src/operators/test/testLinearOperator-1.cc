
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>
#include <cstdlib>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"

#include "materials/Material.h"
#include "operators/LinearOperator.h"
#include "operators/OperatorBuilder.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"

#include "applyTests.h"

void myTest(AMP::UnitTest *ut)
{
  std::string exeName("testLinearOperator-1");
  std::string outerInput_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;
  std::string msgPrefix;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> outerInput_db(new AMP::InputDatabase("outerInput_db"));
  AMP::InputManager::getManager()->parseInputFile(outerInput_file, outerInput_db);
  outerInput_db->printClassData(AMP::plog);

  // Get the Mesh database and create the mesh parameters
  boost::shared_ptr<AMP::Database> database = outerInput_db->getDatabase( "Mesh" );
  boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
  params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

  // Create the meshes from the input database
  AMP::Mesh::Mesh::shared_ptr  meshAdapter = AMP::Mesh::Mesh::buildMesh(params);

  AMP_INSIST( outerInput_db->keyExists("number_of_tests"), "key missing!" );
  int numTests = outerInput_db->getInteger("number_of_tests");

  for(int i = 0; i < numTests; i++) {
    char key[256];
    sprintf(key, "test_%d", i);

    AMP_INSIST( outerInput_db->keyExists(key), "key missing!" );
    std::string innerInput_file = outerInput_db->getString(key);

    boost::shared_ptr<AMP::InputDatabase> innerInput_db(new AMP::InputDatabase("innerInput_db"));
    AMP::InputManager::getManager()->parseInputFile(innerInput_file, innerInput_db);
    innerInput_db->printClassData(AMP::plog);

    AMP_INSIST( innerInput_db->keyExists("testOperator"), "key missing!" );

    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    boost::shared_ptr<AMP::Operator::Operator> testOperator = 
      AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "testOperator",
          innerInput_db,
          elementPhysicsModel);

    msgPrefix=exeName + " : " + innerInput_file;

    if(testOperator.get() != NULL) {
      ut->passes(msgPrefix + " : create");
    } else {
      ut->failure(msgPrefix + " : create");
    }

    boost::shared_ptr<AMP::Operator::LinearOperator> myLinOp =
      boost::dynamic_pointer_cast<AMP::Operator::LinearOperator>(testOperator);

    AMP_INSIST( myLinOp != NULL, "Is not a linear operator!" );

    AMP_INSIST( innerInput_db->keyExists("dofsPerNode"), "key missing!" );
    int dofsPerNode = innerInput_db->getInteger("dofsPerNode");

    AMP::LinearAlgebra::Variable::shared_ptr myInpVar = myLinOp->getInputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr myOutVar = myLinOp->getOutputVariable();
    AMP::Discretization::DOFManager::shared_ptr NodalDOF =
      AMP::Discretization::simpleDOFManager::create(meshAdapter,AMP::Mesh::Vertex,1,dofsPerNode,true);

    {
      AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( NodalDOF, myInpVar, true );
      AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( NodalDOF, myOutVar, true );
      AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( NodalDOF, myOutVar, true );
      // test apply with single variable vectors
      applyTests(ut, msgPrefix, testOperator, rhsVec, solVec, resVec);
    }

    // now run apply tests with multi-vectors
    AMP::LinearAlgebra::Variable::shared_ptr auxInpVar(new AMP::LinearAlgebra::Variable("testLinearOperator-1-auxInpVar"+i));
    AMP::LinearAlgebra::Variable::shared_ptr auxOutVar(new AMP::LinearAlgebra::Variable("testLinearOperator-1-auxOutVar"+i));

    boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> myMultiInpVar( new
        AMP::LinearAlgebra::MultiVariable("MultiInputVariable"));
    myMultiInpVar->add(myInpVar);
    myMultiInpVar->add(auxInpVar);

    boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> myMultiOutVar( new 
        AMP::LinearAlgebra::MultiVariable("MultiOutputVariable"));
    myMultiOutVar->add(myOutVar);
    myMultiOutVar->add(auxOutVar);

    {
      AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( NodalDOF, myMultiInpVar, true );
      AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( NodalDOF, myMultiOutVar, true );
      AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( NodalDOF, myMultiOutVar, true );

      // test apply with single multivariable vectors
      applyTests(ut, msgPrefix, testOperator, rhsVec, solVec, resVec);
    }

    // test getJacobianParameters
    msgPrefix=exeName + " : " + innerInput_file;
    boost::shared_ptr<AMP::LinearAlgebra::Vector> nullGuess;
    boost::shared_ptr<AMP::Operator::OperatorParameters> jacobianParameters = testOperator->getJacobianParameters(nullGuess);

    if(jacobianParameters.get() == NULL) {
      ut->passes(msgPrefix + "getJacobianParameters (should return NULL for now)");
    } else {
      ut->failure(msgPrefix + "getJacobianParameters (should return NULL for now)");
    }

  }//end for i

  ut->passes("testLinearOperator-1");

}

int main(int argc, char *argv[])
{
  AMP::AMPManagerProperties startup_properties;
  startup_properties.use_MPI_Abort = false;
  AMP::AMPManager::startup(argc,argv,startup_properties);
  AMP::UnitTest ut;

  try {
    myTest(&ut);
  } catch (std::exception &err) {
    std::stringstream err_message;
    err_message << "ERROR: Caught standard expection in "<< argv[0] << ": " << err.what();
    std::cout << err_message.str() << std::endl;
    ut.failure(err_message.str());
  } catch( ... ) {
    std::stringstream err_message;
    err_message << "ERROR: Caught expection in "<< argv[0];
    std::cout << err_message.str() << std::endl;
    ut.failure(err_message.str());
  }

  ut.report();

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;
}   


