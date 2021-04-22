// This program simulates the load balance with a given input file on a given number of processors

#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/ampmesh/loadBalance/loadBalanceSimulator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"

#include "ProfilerApp.h"

#include <cmath>
#include <iomanip>
#include <iostream>


// Main function
int run( int N_procs, const std::string &filename, double ratio )
{
    // Simulate loading the mesh
    auto input_db = AMP::Database::parseInputFile( filename );
    auto database = input_db->getDatabase( "Mesh" );
    double t0     = AMP::AMP_MPI::time();
    AMP::Mesh::loadBalanceSimulator mesh( database );
    mesh.setProcs( N_procs );
    double t1 = AMP::AMP_MPI::time();

    // Print the results of the load balance
    mesh.print();

    // Print the time required
    std::cout << "Time = " << t1 - t0 << std::endl << std::endl;

    // Get the worst and average element count
    auto cost  = mesh.getRankCost();
    double max = 0, avg = 0;
    for ( auto x : cost ) {
        avg += x;
        max = std::max( max, x );
    }
    avg /= cost.size();

    // Print the errors and return
    int N_errors = 0;
    if ( t1 - t0 > 12 ) {
        N_errors++;
        std::cout << "load balance failed run time limits" << std::endl;
    }
    if ( max > ratio * avg ) {
        N_errors++;
        std::cout << "load balance failed quality limits" << std::endl;
    }
    return N_errors;
}


// Main function
int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    PROFILE_ENABLE( 3 );
    PROFILE_START( "Main" );

    // Load the inputs
    if ( argc < 3 ) {
        std::cout << "Error calling test_LoadBalancer, format should be:" << std::endl;
        std::cout << "   ./test_LoadBalancer  N_procs  input_file" << std::endl;
        return -1;
    }
    int N_procs   = std::atoi( argv[1] );
    auto filename = argv[2];
    double ratio  = 2.0;
    if ( argc > 3 )
        ratio = atof( argv[3] );

    // Run the problem
    int N_errors = run( N_procs, filename, ratio );

    // Shutdown AMP
    PROFILE_SAVE( filename );
    AMP::AMPManager::shutdown();
    return N_errors;
}
