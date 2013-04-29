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
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "boost/shared_ptr.hpp"


/************************************************************************
*                                                                       *
* This tests whether we can create and use an InputManager object       *
*                                                                       *
************************************************************************/
void mytest(AMP::UnitTest *ut)
{
    std::string input_file = "input_InputManager-1";
    std::string log_file = "output_InputManager-1";

    // Process command line arguments and dump to log file.
    AMP::PIO::logOnlyNodeZero(log_file);

    // Create input database and parse all data in input file.
    boost::shared_ptr<AMP::InputDatabase> input_db ( new AMP::InputDatabase("input_db") );
    AMP::InputManager::getManager()->parseInputFile(input_file, input_db);

    if ( !input_db->getBool("filename") ) 
        ut->failure("InputManager-1");

    input_db.reset();
    
    ut->passes("InputManager read successfully.");
}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    
    try {
        mytest(&ut);
    } catch (std::exception &err) {
        std::cout << "ERROR: While testing test_InputManager-1, " << err.what() << std::endl;
        ut.failure("test_InputManager-1");
    } catch( ... ) {
        std::cout << "ERROR: While testing test_InputManager-1, An unknown exception was thrown." << std::endl;
        ut.failure("test_InputManager-1");
    }
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    if ( num_failed!=0 ) {
        std::cout << "Exiting with errors\n";
        return 1;
    }
    return 0;
}   

