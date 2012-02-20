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

#include "materials/Material.h"
#include "operators/LinearOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"


void myTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST( input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!" );
  std::string mesh_file = input_db->getString("Mesh");

  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter =
    AMP::Mesh::MeshManager::Adapter::shared_ptr(new AMP::Mesh::MeshManager::Adapter());
  meshAdapter->readExodusIIFile ( mesh_file.c_str() );

  AMP_INSIST( input_db->keyExists("testNonlinearMechanicsOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> testNonlinOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
															    "testNonlinearMechanicsOperator",
															    input_db,
															    elementPhysicsModel));
  testNonlinOperator->init();

  AMP::LinearAlgebra::Variable::shared_ptr solVar = testNonlinOperator->getInputVariable(AMP::Operator::Mechanics::DISPLACEMENT); 
  AMP::LinearAlgebra::Variable::shared_ptr resVar = testNonlinOperator->getOutputVariable();

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  AMP::LinearAlgebra::Vector::shared_ptr solVec = meshAdapter->createVector( solVar );
  //AMP::LinearAlgebra::Vector::shared_ptr rhsVec = meshAdapter->createVector( resVar );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = meshAdapter->createVector( resVar );

  solVec->setToScalar(5.0);

  AMP::pout<<"Solution Norm: "<<(solVec->L2Norm())<<std::endl;

  testNonlinOperator->apply(nullVec, solVec, resVec, 1.0, 0.0);

  AMP::pout<<"resNorm1 = "<<(resVec->L2Norm())<<std::endl;

  testNonlinOperator->apply(nullVec, solVec, resVec, 1.0, 0.0);

  AMP::pout<<"resNorm2 = "<<(resVec->L2Norm())<<std::endl;

  ut->passes(exeName);

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back("testNonlinearMechanics-apply-1");

    for(unsigned int i = 0; i < exeNames.size(); i++) {
        try {
            myTest(&ut, exeNames[i]);
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


