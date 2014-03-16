#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "libmesh.h"

#include "operators/ElementOperationFactory.h"

#include "operators/mechanics/MechanicsLinearElement.h"
#include "operators/mechanics/MechanicsNonlinearElement.h"

#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionNonlinearElement.h"

#include "operators/libmesh/MassLinearElement.h"



void ElementOperationFactoryTest(AMP::UnitTest *ut)
{
  std::string exeName("testElementOperationFactory-1");
  std::string outerInput_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> outerInput_db(new AMP::InputDatabase("outerInput_db"));
  AMP::InputManager::getManager()->parseInputFile(outerInput_file, outerInput_db);
  outerInput_db->printClassData(AMP::plog);
  
  AMP_INSIST( outerInput_db->keyExists("number_of_tests"), "key missing!" );
  int numTests = outerInput_db->getInteger("number_of_tests");
  
  for(int i = 0; i < numTests; i++) {
    char key[256];
    sprintf(key, "test_%d", i);
    
    AMP_INSIST( outerInput_db->keyExists(key), "key missing!" );
    std::string inputFile = outerInput_db->getString(key);
    boost::shared_ptr<AMP::InputDatabase> innerInput_db(new AMP::InputDatabase("innerInput_db"));
    AMP::InputManager::getManager()->parseInputFile(inputFile, innerInput_db);
    innerInput_db->printClassData(AMP::plog);
    
    AMP_INSIST(innerInput_db->keyExists("ElementOperation"), "Key ''ElementOperation'' is missing!");
    boost::shared_ptr<AMP::Database> elemOp_db = innerInput_db->getDatabase("ElementOperation");
    std::string name = elemOp_db->getString("name");
    boost::shared_ptr<AMP::Operator::ElementOperation> elementOperation =
      AMP::Operator::ElementOperationFactory::createElementOperation(elemOp_db);
    
    if(name=="MechanicsLinearElement")
      {
    boost::shared_ptr<AMP::Operator::MechanicsLinearElement> mechOperation = boost::dynamic_pointer_cast<AMP::Operator::MechanicsLinearElement>(elementOperation );
    if(mechOperation.get()!=NULL)
      {
        ut->passes(exeName + " : " + inputFile + " : MechanicsLinearElement");
      }
    else
      {
        ut->failure(exeName + " : " + inputFile + " : MechanicsLinearElement");
      }
      }
    else if(name=="DiffusionLinearElement")
      {
    boost::shared_ptr<AMP::Operator::DiffusionLinearElement> diffusionOperation = boost::dynamic_pointer_cast<AMP::Operator::DiffusionLinearElement>(elementOperation);
    if(diffusionOperation.get()!=NULL)
      {
        ut->passes(exeName + " : " + inputFile + " : DiffusionLinearElement");
        }
    else
      {
        ut->failure(exeName + " : " + inputFile + " : DiffusionLinearElement");
      }
      }
    else if(name=="MechanicsNonlinearElement")
      {
    boost::shared_ptr<AMP::Operator::MechanicsNonlinearElement> mechOperation = boost::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearElement>(elementOperation );
    if(mechOperation.get()!=NULL)
      {
        ut->passes(exeName + " : " + inputFile + " : MechanicsNonlinearElement");
      }
    else
      {
        ut->failure(exeName + " : " + inputFile + " : MechanicsNonlinearElement");
      }
      }
    else if(name=="DiffusionNonlinearElement")
      {
    boost::shared_ptr<AMP::Operator::DiffusionNonlinearElement> diffusionOperation = boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearElement>(elementOperation);
    if(diffusionOperation.get()!=NULL)
      {
        ut->passes(exeName + " : " + inputFile + " : DiffusionNonlinearElement");
      }
    else
      {
        ut->failure(exeName + " : " + inputFile + " : DiffusionNonlinearElement");
      }
      }
    else if(name=="MassLinearElement")
      {
    boost::shared_ptr<AMP::Operator::MassLinearElement> massOperation = boost::dynamic_pointer_cast<AMP::Operator::MassLinearElement>(elementOperation);
    if(massOperation.get()!=NULL)
      {
        ut->passes(exeName + " : " + inputFile + " : MassLinearElement");
        }
    else
        {
          ut->failure(exeName + " : " + inputFile + " : MassLinearElement");
        }
      }
  }

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
 
    try {
        ElementOperationFactoryTest(&ut);
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


