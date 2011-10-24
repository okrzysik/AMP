#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include <string>

#include <assert.h>

#include <fstream>

#include <sys/stat.h>
#include "../Database.h"
#include "../InputDatabase.h"
#include "../InputManager.h"
#include "../AMP_MPI.h"
#include "../AMPManager.h"
#include "../PIO.h"
#include "boost/shared_ptr.hpp"


/************************************************************************
*                                                                       *
* This tests whether we can create and use an InputManager object       *
*                                                                       *
************************************************************************/
void mytest(AMP::UnitTest *ut)
{
    std::string input_file = "input_Database-1";
    std::string log_file = "output_Database-1";

    // Process command line arguments and dump to log file.
    std::string info = "test_Database-1: tests the getDatabase function and the destructor of AMP::Database";
    AMP::PIO::logOnlyNodeZero(log_file);

    // Create input database and parse all data in input file.
    boost::shared_ptr<AMP::InputDatabase> input_db ( new AMP::InputDatabase("input_db") );
    AMP::InputManager::getManager()->parseInputFile(input_file, input_db);

    boost::shared_ptr<AMP::Database> tmp_db = input_db->getDatabase("Try");
    int number = tmp_db->getInteger("number");

    if(number>0) {
      std::vector<int> intArray = tmp_db->getIntegerArray("intArray");
      if( (int) intArray.size() != number ) ut->failure("intArray was the wrong size");
      std::vector<double> doubleArray = tmp_db->getDoubleArray("doubleArray");
      if( (int) doubleArray.size() != number ) ut->failure("doubleArray was the wrong size");
    }

    std::cout<<"tmp_db created and destroyed successfully."<<std::endl;

    input_db.reset();
    
    ut->passes("Get Database Works and the Destructor of AMP::Database works.");
}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    
    try {
        mytest(&ut);
    } catch (std::exception &err) {
        std::cout << "ERROR: While testing test_Database-1, " << err.what() << std::endl;
        ut.failure("test_Database-1");
    } catch( ... ) {
        std::cout << "ERROR: While testing test_InputManager-1, An unknown exception was thrown." << std::endl;
        ut.failure("test_Database-1");
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   

