
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

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshVariable.h"

#include "libmesh.h"

#include "materials/Material.h"
#include "../LinearOperator.h"
#include "../ColumnOperator.h"
#include "../OperatorBuilder.h"

#include "applyTests.h"

void myTest(AMP::UnitTest *ut)
{
  std::string exeName("testLinearColumnOperator-1");
  std::string outerInput_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;
  std::string msgPrefix;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> outerInput_db(new AMP::InputDatabase(
        "outerInput_db"));
  AMP::InputManager::getManager()->parseInputFile(outerInput_file,
      outerInput_db);
  outerInput_db->printClassData(AMP::plog);

  AMP_INSIST(outerInput_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  std::string mesh_file = outerInput_db->getString("Mesh");

  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter =
    AMP::Mesh::MeshManager::Adapter::shared_ptr(
        new AMP::Mesh::MeshManager::Adapter());
  meshAdapter->readExodusIIFile(mesh_file.c_str());

  AMP_INSIST(outerInput_db->keyExists("number_of_tests"), "key missing!");
  int numTests = outerInput_db->getInteger("number_of_tests");

  for (int i = 0; i < numTests; i++)
  {
    char key[256];
    sprintf(key, "test_%d", i);

    AMP_INSIST(outerInput_db->keyExists(key), "key missing!");
    std::string innerInput_file = outerInput_db->getString(key);

    boost::shared_ptr<AMP::InputDatabase> innerInput_db(
        new AMP::InputDatabase("innerInput_db"));
    AMP::InputManager::getManager()->parseInputFile(innerInput_file,
        innerInput_db);
    innerInput_db->printClassData(AMP::plog);

    AMP_INSIST(innerInput_db->keyExists("numberOfOperators"), "key missing!");

    const int numberOfOperators = innerInput_db->getInteger(
        "numberOfOperators");

    AMP_INSIST(numberOfOperators >= 1,
        "more than zero operators need to be specified");

    // create a column operator object
    boost::shared_ptr<AMP::Operator::OperatorParameters> params;
    boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator(
        new AMP::Operator::ColumnOperator(params));

    boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> columnInputVariable(
        new AMP::LinearAlgebra::MultiVariable("columnInputVariable"));

    double defTemp = -1.0;
    double defConc = -1.0;
    size_t nVars=0;
    for (int opN = 1; opN <= numberOfOperators; opN++)
    {
      char testOpName[256];
      sprintf(testOpName, "testOperator%d", opN);
      AMP_INSIST(innerInput_db->keyExists(testOpName), "key missing!");

      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
      boost::shared_ptr<AMP::Database> testOp_db =
        innerInput_db->getDatabase(testOpName);
      boost::shared_ptr<AMP::Operator::Operator> testOperator =
        AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
						       testOpName,
						       innerInput_db,
						       elementPhysicsModel);

      boost::shared_ptr<AMP::Operator::LinearOperator> myLinOp =
        boost::dynamic_pointer_cast<AMP::Operator::LinearOperator>(
            testOperator);
      AMP_INSIST(myLinOp != NULL, "Is not a linear operator!");

      columnOperator->append(testOperator);

      AMP::LinearAlgebra::Variable::shared_ptr opVar = myLinOp->getInputVariable();

      if (opVar.get() != NULL)
      {
        columnInputVariable->add(opVar);
        nVars++;
      }

      // this only works as long at least one of the operators is diffusion and
      // its transport model has defaults defined
      boost::shared_ptr<AMP::Database> model_db;
      if (testOp_db->keyExists("VolumeOperator"))
	{
	  boost::shared_ptr<AMP::Database> volOp_db = innerInput_db->getDatabase(testOp_db->getString("VolumeOperator"));
	  if ((volOp_db->getName() == "DiffusionNonlinearFEOperator")||(volOp_db->getName() == "DiffusionLinearFEOperator"))
	    {
	      model_db = innerInput_db->getDatabase(volOp_db->getString("LocalModel"));
	    }
        }

      if (model_db) {
        defTemp = model_db->getDouble("Default_Temperature");
        defConc = model_db->getDouble("Default_Concentration");
      }
    }

    msgPrefix = exeName + " : " + innerInput_file;
    ut->passes(msgPrefix + " : create");

    boost::shared_ptr<AMP::Operator::Operator> testOperator =
      boost::dynamic_pointer_cast<AMP::Operator::Operator>(columnOperator);

    AMP::LinearAlgebra::Variable::shared_ptr columnOutputVariable =
      columnOperator->getOutputVariable();

    {
      AMP::LinearAlgebra::Vector::shared_ptr solVec = meshAdapter->createVector(
          columnInputVariable);
      AMP::LinearAlgebra::Vector::shared_ptr rhsVec = meshAdapter->createVector(
          columnOutputVariable);
      AMP::LinearAlgebra::Vector::shared_ptr resVec = meshAdapter->createVector(
          columnOutputVariable);

      for (size_t i=0; i<nVars; i++) {
        AMP::LinearAlgebra::Variable::shared_ptr opVar = columnInputVariable->getVariable(i);
        if (opVar->getName()=="temperature") {
          AMP::LinearAlgebra::Vector::shared_ptr tVec = solVec->subsetVectorForVariable(columnInputVariable->getVariable(i));
          tVec->setToScalar(defTemp);
        }
        if (opVar->getName()=="concentration") {
          AMP::LinearAlgebra::Vector::shared_ptr cVec = solVec->subsetVectorForVariable(columnInputVariable->getVariable(i));
          cVec->setToScalar(defConc);
        }
      }

      // test apply with single variable vectors
      applyTests(ut, msgPrefix, testOperator, rhsVec, solVec, resVec);

    }
#if 0
    // test getJacobianParameters
    msgPrefix=exeName + " : " + innerInput_file;
    boost::shared_ptr<AMP::LinearAlgebra::Vector> nullGuess;
    boost::shared_ptr<AMP::Operator::OperatorParameters> jacobianParameters = testOperator->getJacobianParameters(nullGuess);
    if(jacobianParameters.get()!=NULL)
    {
      ut.passes(msgPrefix + "getJacobianParameters (should return NULL for now)");
    }
    else
    {
      ut.numFails++;
    }
#endif

  }//end for i

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

