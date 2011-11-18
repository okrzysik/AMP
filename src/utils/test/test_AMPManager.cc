#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"


//  This test will start and shutdown AMP
int main(int argc, char *argv[])
{
    #ifdef USE_MPI
        // Initialize MPI directly
        MPI_Init(&argc, &argv);
    #endif
    // Control the behavior of the startup
    AMP::AMPManagerProperties startup_properties;
    startup_properties.print_times = true;
    startup_properties.use_MPI_Abort = false;
    // Start AMP
    AMP::AMPManager::startup(argc,argv,startup_properties);
    int globalSize=0;
    // Limit the scope of variables
    { 
        // Create the global comm
        AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
        globalSize = globalComm.getSize();
        // Introduce a memory leak to catch in valgrind later
        double *x = new double[100];
        if ( x==NULL ) 
            AMP_ERROR("error");
        // Test the abort
        try {
            AMP_ERROR("Catch this");
            return -1;
        } catch (...) {
            // This is correct
        }
    }
    // Shutdown
    AMP::AMPManager::shutdown();
    // If we are using MPI, try to initailize AMPManager on a different communicator
    #ifdef USE_MPI
        int size, rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Comm_size(MPI_COMM_WORLD,&size);                
        if ( size!=globalSize )
            return -1;
        // Create a new comm to initialize AMPManager
        MPI_Comm new_comm_world;
        int color = rank%2;
        MPI_Comm_split(MPI_COMM_WORLD,color,0,&new_comm_world);
        // Start AMP
        startup_properties.COMM_WORLD = new_comm_world;
        AMP::AMPManager::startup(argc,argv,startup_properties);
        // Limit the scope of variables
        { 
            AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
            int newGlobalSize = globalComm.getSize();
            if ( newGlobalSize==globalSize && globalSize!=1 )
                return -1;
        }
        // Shutdown
        AMP::AMPManager::shutdown();
    #endif
    #ifdef USE_MPI
        // Close MPI
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    #endif
    // Finished successfully
    return 0;
}   

