// This file tests the HDF5 interfaces
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <random>
#include <string>

#include "AMP/IO/HDF5.h"
#include "AMP/IO/HDF5_Class.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"


// Main
int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );

    if ( argc == 1 ) {
        std::cerr << "view_HDF5 filename" << std::endl;
        return -1;
    }

    // Loop through the input files loading and printing the variables
    for ( int i = 1; i < argc; i++ ) {
        std::cout << argv[i] << ":\n";
        if ( !AMP::Utilities::fileExists( argv[i] ) ) {
            std::cerr << "File does not exist\n";
            return -1;
        }
        auto fid  = AMP::openHDF5( argv[i], "r" );
        auto data = AMP::readHDF5( fid, "/" );
        AMP::closeHDF5( fid );
        data->print( 2, "  " );
        std::cout << std::endl;
    }

    AMP::AMPManager::shutdown();
    return 0;
}
