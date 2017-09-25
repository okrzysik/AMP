// This program simulates the load balance with a given input file on a given number of processors

#include "ampmesh/Mesh.h"
#include "ampmesh/loadBalance/loadBalanceSimulator.h"
#include "utils/AMPManager.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"

#include <iomanip>
#include <iostream>
#include <math.h>


// Main function
int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );

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
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    AMP::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    std::vector<int> comm_ranks( N_procs );
    for ( int i = 0; i < N_procs; i++ )
        comm_ranks[i] = i;
    double t0 = AMP::AMP_MPI::time();
    AMP::Mesh::loadBalanceSimulator mesh( params, comm_ranks );
    double t1 = AMP::AMP_MPI::time();

    // Print the results of the load balance
    if ( N_procs < 10000 ) {
        mesh.print();
        std::cout << std::endl;
    }
    std::cout << "min = " << mesh.min() << std::endl;
    std::cout << "max = " << mesh.max() << std::endl;
    std::cout << "avg = " << mesh.avg() << std::endl;
    std::cout << "time = " << t1 - t0 << std::endl;

    // Shutdown AMP
    AMP::AMPManager::shutdown();

    // Print the errors and return
    int N_errors = 0;
    if ( t1 - t0 > 10 ) {
        N_errors++;
        std::cout << "load balance failed run time limits" << std::endl;
    }
    if ( ( (double) mesh.max() ) > ratio * ( (double) mesh.avg() ) ) {
        N_errors++;
        std::cout << "load balance failed quality limits" << std::endl;
    }
    return N_errors;
}
