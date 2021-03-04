// This file tests the HDF5 interfaces
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <random>
#include <string>

#include "AMP/utils/AMPManager.h"
#include "AMP/utils/HDF5_Class.h"
#include "AMP/utils/HDF5_IO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"


// Main
int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );

    if ( argc != 2 ) {
        std::cerr << "view_HDF5 filename" << std::endl;
        return -1;
    }

    // Open the HDF5 file
    auto fid = AMP::openHDF5( argv[1], "r" );

    // Read the data
    auto data = AMP::readHDF5( fid, "/" );

    // Print the data
    data->print( 2, "  " );

    data.reset();
    AMP::AMPManager::shutdown();
    return 0;
}
