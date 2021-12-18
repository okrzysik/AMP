#include "AMP/IO/PIO.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include <memory>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>


void mytest( AMP::UnitTest *ut )
{
    std::string input_file = "input_InputManager";
    std::string log_file   = "output_InputManager";

    // Create input database and parse all data in input file
    auto input_db = AMP::Database::parseInputFile( input_file );

    if ( !input_db->getScalar<bool>( "filename" ) )
        ut->failure( "InputManager-1" );

    input_db.reset();

    ut->passes( "InputManager read successfully." );
}


int main( int argc, char *argv[] )
{
    AMP::AMP_MPI::start_MPI( argc, argv );
    AMP::UnitTest ut;

    mytest( &ut );
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMP_MPI::stop_MPI();
    if ( num_failed != 0 ) {
        std::cout << "Exiting with errors\n";
        return 1;
    }
    return 0;
}
