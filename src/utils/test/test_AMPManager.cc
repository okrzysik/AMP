#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "utils/AMPManager.h"


//  This test will start and shutdown AMP
int main(int argc, char *argv[])
{
    // Control the behavior of the startup
    AMP::AMPManagerProperties startup_properties;
    startup_properties.print_times = true;
    // Start AMP
    AMP::AMPManager::startup(argc,argv,startup_properties);
    // Introduce a memory leak to catch later
    double *x = new double[100];
    if ( x==NULL ) 
        AMP_ERROR("error");
    // Shutdown
    AMP::AMPManager::shutdown();
    return 0;
}   

