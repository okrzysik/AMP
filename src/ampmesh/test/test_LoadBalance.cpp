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


// Helper functions
static inline size_t min( const std::vector<double> &x )
{
    double y = 1e100;
    for ( auto i : x )
        y = std::min( y, i );
    return y;
}
static inline size_t max( const std::vector<double> &x )
{
    double y = 0;
    for ( auto i : x )
        y = std::max( y, i );
    return y;
}
static inline size_t avg( const std::vector<double> &x )
{
    double y = 0;
    for ( auto i : x )
        y += i;
    return y / x.size();
}


// Main function
int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    PROFILE_ENABLE( 3 );
    PROFILE_START( "Main" );

    // Load the input file and the desired number of processors
    if ( argc < 3 ) {
        std::cout << "Error calling test_LoadBalancer, format should be:" << std::endl;
        std::cout << "   ./test_LoadBalancer  N_procs  input_file" << std::endl;
        return -1;
    }
    int N_procs            = std::atoi( argv[1] );
    std::string input_file = argv[2];
    double ratio           = 2.0;
    if ( argc > 3 )
        ratio = atof( argv[3] );

    // Simulate loading the mesh
    auto input_db = AMP::Database::parseInputFile( input_file );
    auto database = input_db->getDatabase( "Mesh" );
    double t0     = AMP::AMP_MPI::time();
    AMP::Mesh::loadBalanceSimulator mesh( database );
    mesh.setProcs( N_procs );
    double t1 = AMP::AMP_MPI::time();

    // Print the results of the load balance
    if ( N_procs < 10000 ) {
        mesh.print();
        std::cout << std::endl;
    }
    auto cost = mesh.getRankCost();
    std::cout << "min = " << min( cost ) << std::endl;
    std::cout << "max = " << max( cost ) << std::endl;
    std::cout << "avg = " << avg( cost ) << std::endl;
    std::cout << "time = " << t1 - t0 << std::endl;

    // Shutdown AMP
    PROFILE_SAVE( input_file );
    AMP::AMPManager::shutdown();

    // Print the errors and return
    int N_errors = 0;
    if ( t1 - t0 > 10 ) {
        N_errors++;
        std::cout << "load balance failed run time limits" << std::endl;
    }
    if ( max( cost ) > ratio * avg( cost ) ) {
        N_errors++;
        std::cout << "load balance failed quality limits" << std::endl;
    }
    return N_errors;
}
